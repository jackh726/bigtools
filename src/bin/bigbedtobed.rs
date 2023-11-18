use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::Path;

use bigtools::utils::streaming_linereader::StreamingLineReader;
use clap::Parser;

use futures::FutureExt;

use bigtools::utils::reopen::{Reopen, SeekableRead};
use bigtools::utils::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};
use bigtools::{BBIReadError, BigBedRead, ChromInfo};
use tokio::runtime;

pub fn write_bed<R: Reopen + SeekableRead + Send + 'static>(
    bigbed: BigBedRead<R>,
    mut out_file: File,
    nthreads: usize,
    chrom: Option<String>,
    start: Option<u32>,
    end: Option<u32>,
) -> Result<(), BBIReadError> {
    let runtime = if nthreads == 1 {
        runtime::Builder::new_current_thread().build().unwrap()
    } else {
        runtime::Builder::new_multi_thread()
            .worker_threads(nthreads)
            .build()
            .unwrap()
    };

    let chrom_files: Vec<io::Result<(_, TempFileBuffer<File>)>> = bigbed
        .chroms()
        .into_iter()
        .cloned()
        .filter(|c| chrom.as_ref().map_or(true, |chrom| &c.name == chrom))
        .map(|chrom| {
            let bigbed = bigbed.reopen()?;
            let (buf, file): (TempFileBuffer<File>, TempFileBufferWriter<File>) =
                TempFileBuffer::new();
            let writer = io::BufWriter::new(file);
            async fn file_future<R: SeekableRead + 'static>(
                mut bigbed: BigBedRead<R>,
                chrom: ChromInfo,
                mut writer: io::BufWriter<TempFileBufferWriter<File>>,
                start: Option<u32>,
                end: Option<u32>,
            ) -> Result<(), BBIReadError> {
                for raw_val in bigbed.get_interval(&chrom.name, 0, chrom.length)? {
                    let val = raw_val?;
                    if let Some(start) = start {
                        if val.start <= start {
                            continue;
                        }
                    }
                    if let Some(end) = end {
                        if val.start > end {
                            continue;
                        }
                    }
                    let end = if !val.rest.is_empty() {
                        format!("\t{}\n", val.rest)
                    } else {
                        "\n".to_string()
                    };
                    writer.write_fmt(format_args!(
                        "{}\t{}\t{}{}",
                        chrom.name, val.start, val.end, end
                    ))?;
                }
                Ok(())
            }
            let handle = runtime
                .spawn(file_future(bigbed, chrom, writer, start, end))
                .map(|f| f.unwrap());
            Ok((handle, buf))
        })
        .collect::<Vec<_>>();

    for res in chrom_files {
        let (f, mut buf) = res.unwrap();
        buf.switch(out_file);
        runtime.block_on(f).unwrap();
        out_file = buf.await_real_file();
    }

    Ok(())
}

pub fn write_bed_from_bed<R: Reopen + SeekableRead + Send + 'static>(
    mut bigbed: BigBedRead<R>,
    out_file: File,
    bed: File,
) -> Result<(), BBIReadError> {
    let mut bedstream = StreamingLineReader::new(BufReader::new(bed));
    let mut writer = io::BufWriter::new(out_file);

    while let Some(line) = bedstream.read() {
        let line = line?;
        let mut split = line.trim().splitn(5, '\t');
        let chrom = split.next().expect("Missing chrom");
        let start = split.next().expect("Missing start").parse::<u32>().unwrap();
        let end = split.next().expect("Missing end").parse::<u32>().unwrap();
        for raw_val in bigbed.get_interval(chrom, start, end)? {
            let val = raw_val?;
            let end = if !val.rest.is_empty() {
                format!("\t{}\n", val.rest)
            } else {
                "\n".to_string()
            };
            writer.write_fmt(format_args!("{}\t{}\t{}{}", chrom, val.start, val.end, end))?;
        }
    }

    Ok(())
}

#[derive(Parser)]
#[command(about = "Converts an input bigBed to a bed. Can be multi-threaded for substantial speedups. Note for roughly each core, one temporary file will be opened.", long_about = None)]
struct Cli {
    /// the bigbed to get convert to bed
    big_bed: String,

    /// the path of the bed to output to
    bed: String,

    /// If set, restrict output to given chromosome
    #[arg(long)]
    chrom: Option<String>,

    /// If set, restrict output to regions greater than or equal to it
    #[arg(long)]
    start: Option<u32>,

    /// If set, restrict output to regions less than it
    #[arg(long)]
    end: Option<u32>,

    /// If set, restrict output to regions overlapping the bed file
    overlap_bed: Option<String>,

    /// Set the number of threads to use. This tool will nearly always benefit from more cores (<= # chroms). Note: for parts of the runtime, the actual usage may be nthreads+1
    #[arg(short = 't', long)]
    #[arg(default_value_t = 6)]
    nthreads: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args_os().map(|a| {
        bigtools::compat_replace!(a;
            replace:
                "-chrom", "--chrom";
                "-start", "--start";
                "-end", "--end";
                "-bed", "-overlap-bed"
            ignore:
            unimplemented:
                "-header";
                "-udcDir"
        )
    });
    let matches = Cli::parse_from(args);

    let bigbedpath = matches.big_bed;
    let bedpath = matches.bed;

    let nthreads = matches.nthreads;

    let bigbed = BigBedRead::open_file(&bigbedpath)?;
    let bed = File::create(bedpath)?;

    if matches.start.is_some() || matches.end.is_some() & matches.chrom.is_none() {
        eprintln!("Cannot specify --start or --end without specifying --chrom.");
        return Ok(());
    }

    if matches.chrom.is_some() && matches.overlap_bed.is_some() {
        eprintln!("Cannot specify both --overlap-bed and interval to overlap.");
        return Ok(());
    }

    match matches.overlap_bed {
        Some(overlap_bed) => {
            if !Path::exists(&Path::new(&overlap_bed)) {
                eprintln!("Overlap bed file does not exist.");
                return Ok(());
            }
            let overlap_bed = File::open(overlap_bed)?;
            write_bed_from_bed(bigbed, bed, overlap_bed)?;
        }
        None => {
            write_bed(
                bigbed,
                bed,
                nthreads,
                matches.chrom,
                matches.start,
                matches.end,
            )?;
        }
    }

    Ok(())
}

#[test]
fn verify_cli_bigbedtobed() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
