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
use bigtools::{BBIReadError, BigWigRead, ChromInfo};
use tokio::runtime;
use ufmt::uwrite;

pub fn write_bg_singlethreaded<R: SeekableRead + Send + 'static>(
    mut bigwig: BigWigRead<R>,
    out_file: File,
    chrom: Option<String>,
    start: Option<u32>,
    end: Option<u32>,
) -> Result<(), BBIReadError> {
    let start = chrom.as_ref().and_then(|_| start);
    let end = chrom.as_ref().and_then(|_| end);

    let mut chroms: Vec<ChromInfo> = bigwig.chroms().to_vec();
    chroms.sort_by(|a, b| a.name.cmp(&b.name));
    let mut writer = io::BufWriter::with_capacity(32 * 1000, out_file);
    for chrom in chroms {
        let start = start.unwrap_or(0);
        let end = end.unwrap_or(chrom.length);
        let mut values = bigwig.get_interval(&chrom.name, start, end)?;
        while let Some(raw_val) = values.next() {
            let val = raw_val?;

            let mut buf = String::with_capacity(50); // Estimate
            uwrite!(
                &mut buf,
                "{}\t{}\t{}\t{}\n",
                chrom.name,
                val.start,
                val.end,
                ryu::Buffer::new().format(val.value)
            )
            .unwrap();
            writer.write(buf.as_bytes())?;
        }
    }

    Ok(())
}

pub async fn write_bg<R: Reopen + SeekableRead + Send + 'static>(
    bigwig: BigWigRead<R>,
    mut out_file: File,
    chrom: Option<String>,
    start: Option<u32>,
    end: Option<u32>,
    runtime: &runtime::Handle,
) -> Result<(), BBIReadError> {
    let start = chrom.as_ref().and_then(|_| start);
    let end = chrom.as_ref().and_then(|_| end);

    let chrom_files: Vec<io::Result<(_, TempFileBuffer<File>)>> = bigwig
        .chroms()
        .into_iter()
        .cloned()
        .filter(|c| chrom.as_ref().map_or(true, |chrom| &c.name == chrom))
        .map(|chrom| {
            let bigwig = bigwig.reopen()?;
            let (buf, file): (TempFileBuffer<File>, TempFileBufferWriter<File>) =
                TempFileBuffer::new();
            let writer = io::BufWriter::new(file);
            async fn file_future<R: Reopen + SeekableRead + 'static>(
                mut bigwig: BigWigRead<R>,
                chrom: ChromInfo,
                mut writer: io::BufWriter<TempFileBufferWriter<File>>,
                start: u32,
                end: u32,
            ) -> Result<(), BBIReadError> {
                for raw_val in bigwig.get_interval(&chrom.name, start, end)? {
                    let val = raw_val?;
                    let mut buf = String::with_capacity(50); // Estimate

                    // Using ryu for f32 to string conversion has a ~15% speedup
                    uwrite!(
                        &mut buf,
                        "{}\t{}\t{}\t{}\n",
                        chrom.name.as_str(),
                        val.start,
                        val.end,
                        ryu::Buffer::new().format(val.value)
                    )
                    .unwrap();
                    writer.write(buf.as_bytes())?;
                }
                Ok(())
            }
            let start = start.unwrap_or(0);
            let end = end.unwrap_or(chrom.length);
            let handle = runtime
                .spawn(file_future(bigwig, chrom, writer, start, end))
                .map(|f| f.unwrap());
            Ok((handle, buf))
        })
        .collect::<Vec<_>>();

    for res in chrom_files {
        let (f, mut buf) = res.unwrap();
        buf.switch(out_file);
        f.await.unwrap();
        while !buf.is_real_file_ready() {
            tokio::task::yield_now().await;
        }
        out_file = buf.await_real_file();
    }

    Ok(())
}

pub fn write_bg_from_bed<R: Reopen + SeekableRead + Send + 'static>(
    mut bigbed: BigWigRead<R>,
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
            let mut buf = String::with_capacity(50); // Estimate

            // Using ryu for f32 to string conversion has a ~15% speedup
            uwrite!(
                &mut buf,
                "{}\t{}\t{}\t{}\n",
                chrom,
                val.start,
                val.end,
                ryu::Buffer::new().format(val.value)
            )
            .unwrap();
            writer.write(buf.as_bytes())?;
        }
    }

    Ok(())
}

#[derive(Parser)]
#[command(about = "Converts an input bigWig to a bedGraph. Can be multi-threaded for substantial speedups. Note for roughly each core, one temporary file will be opened.", long_about = None)]
struct Cli {
    /// the bigwig to get convert to bedgraph
    bigwig: String,

    /// the path of the bedgraph to output to
    bedgraph: String,

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
                "-end", "--end"
            ignore:
            unimplemented:
                "-udcDir"
        )
    });
    let matches = Cli::parse_from(args);

    let bigwigpath = matches.bigwig;
    let bedgraphpath = matches.bedgraph;

    let nthreads = matches.nthreads;

    let bigwig = BigWigRead::open_file(&bigwigpath)?;
    let bedgraph = File::create(bedgraphpath)?;

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
            write_bg_from_bed(bigwig, bedgraph, overlap_bed)?;
        }
        None => {
            if nthreads == 1 {
                write_bg_singlethreaded(
                    bigwig,
                    bedgraph,
                    matches.chrom,
                    matches.start,
                    matches.end,
                )?;
            } else {
                let runtime = runtime::Builder::new_multi_thread()
                    .worker_threads(nthreads)
                    .build()
                    .unwrap();

                runtime.block_on(write_bg(
                    bigwig,
                    bedgraph,
                    matches.chrom,
                    matches.start,
                    matches.end,
                    runtime.handle(),
                ))?;
            }
        }
    }

    dbg!();

    Ok(())
}

#[test]
fn verify_cli_bigwigtobedgraph() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
