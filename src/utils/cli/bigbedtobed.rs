use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::Path;

use clap::Parser;
use futures::FutureExt;
use tokio::runtime;

use crate::utils::reopen::{Reopen, SeekableRead};
use crate::utils::streaming_linereader::StreamingLineReader;
use crate::utils::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};
use crate::{BBIReadError, BigBedRead, ChromInfo};

#[derive(Debug, Parser)]
#[command(
    name = "bigbedtobed",
    about = "Converts an input bigBed to a bed.",
    long_about = "Converts an input bigBed to a bed. Can be multi-threaded for substantial speedups. Note for roughly each core, one temporary file will be opened."
)]
pub struct BigBedToBedArgs {
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

pub fn bigbedtobed(args: BigBedToBedArgs) -> Result<(), Box<dyn Error>> {
    let bigbedpath = args.big_bed;
    let bedpath = args.bed;

    let nthreads = args.nthreads;

    let bigbed = BigBedRead::open_file(&bigbedpath)?;
    let bed = File::create(bedpath)?;

    if args.start.is_some() || args.end.is_some() & args.chrom.is_none() {
        eprintln!("Cannot specify --start or --end without specifying --chrom.");
        return Ok(());
    }

    if args.chrom.is_some() && args.overlap_bed.is_some() {
        eprintln!("Cannot specify both --overlap-bed and interval to overlap.");
        return Ok(());
    }

    match args.overlap_bed {
        Some(overlap_bed) => {
            if !Path::exists(&Path::new(&overlap_bed)) {
                eprintln!("Overlap bed file does not exist.");
                return Ok(());
            }
            let overlap_bed = File::open(overlap_bed)?;
            write_bed_from_bed(bigbed, bed, overlap_bed)?;
        }
        None => {
            write_bed(bigbed, bed, nthreads, args.chrom, args.start, args.end)?;
        }
    }

    Ok(())
}

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
