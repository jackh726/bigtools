use std::collections::VecDeque;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::Path;

use crate::utils::streaming_linereader::StreamingLineReader;
use clap::Parser;

use futures::FutureExt;

use crate::utils::reopen::{Reopen, SeekableRead};
use crate::utils::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};
use crate::{BBIReadError, BigWigRead, ChromInfo};
use tokio::runtime;
use ufmt::uwrite;

#[derive(Clone, Debug, PartialEq, Parser)]
#[command(
    name = "bigwigtobedgraph",
    about = "Converts an input bigWig to a bedGraph.",
    long_about = "Converts an input bigWig to a bedGraph. Can be multi-threaded for substantial speedups."
)]
pub struct BigWigToBedGraphArgs {
    /// the bigwig to get convert to bedgraph
    pub bigwig: String,

    /// the path of the bedgraph to output to
    pub bedgraph: String,

    /// If set, restrict output to given chromosome
    #[arg(long)]
    pub chrom: Option<String>,

    /// If set, restrict output to regions greater than or equal to it
    #[arg(long)]
    pub start: Option<u32>,

    /// If set, restrict output to regions less than it
    #[arg(long)]
    pub end: Option<u32>,

    /// If set, restrict output to regions overlapping the bed file
    pub overlap_bed: Option<String>,

    /// Set the number of threads to use. This tool will nearly always benefit from more cores (<= # chroms).
    #[arg(short = 't', long)]
    #[arg(default_value_t = 6)]
    pub nthreads: usize,

    /// Do not create temporary files for intermediate data. (Only applicable when using multiple threads.)
    /// By default, approximately one temporary file will be opened for each core.
    #[arg(long)]
    #[arg(default_value_t = false)]
    pub inmemory: bool,
}

pub fn bigwigtobedgraph(args: BigWigToBedGraphArgs) -> Result<(), Box<dyn Error>> {
    let bigwigpath = args.bigwig;
    let bedgraphpath = args.bedgraph;

    let nthreads = args.nthreads;

    let bigwig = BigWigRead::open_file(&bigwigpath)?;
    let bedgraph = File::create(bedgraphpath)?;

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
            write_bg_from_bed(bigwig, bedgraph, overlap_bed)?;
        }
        None => {
            // Right now, we don't offload decompression to separate threads,
            // so specifying `chrom` effectively means single-threaded
            if nthreads == 1 || args.chrom.is_some() {
                write_bg_singlethreaded(bigwig, bedgraph, args.chrom, args.start, args.end)?;
            } else {
                write_bg(bigwig, bedgraph, args.inmemory, nthreads)?;
            }
        }
    }

    Ok(())
}

pub fn write_bg_singlethreaded<R: SeekableRead + Send + 'static>(
    mut bigwig: BigWigRead<R>,
    out_file: File,
    chrom: Option<String>,
    start: Option<u32>,
    end: Option<u32>,
) -> Result<(), BBIReadError> {
    let start = chrom.as_ref().and_then(|_| start);
    let end = chrom.as_ref().and_then(|_| end);

    let chroms: Vec<ChromInfo> = if let Some(arg_chrom) = chrom {
        let chrom = bigwig.chroms().iter().find(|c| c.name == arg_chrom);
        let Some(chrom) = chrom else {
            eprintln!("{arg_chrom} not found in file.");
            return Ok(());
        };
        vec![chrom.clone()]
    } else {
        bigwig.chroms().to_vec()
    };
    let mut writer = io::BufWriter::with_capacity(32 * 1000, out_file);
    for chrom in chroms {
        let start = start.unwrap_or(0);
        let end = end.unwrap_or(chrom.length);
        let mut values = bigwig.get_interval(&chrom.name, start, end)?;
        let mut buf = String::with_capacity(50); // Estimate
        while let Some(raw_val) = values.next() {
            let val = raw_val?;

            // Using ryu for f32 to string conversion has a ~15% speedup
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
            buf.clear();
        }
    }

    Ok(())
}

pub fn write_bg<R: Reopen + SeekableRead + Send + 'static>(
    bigwig: BigWigRead<R>,
    mut out_file: File,
    inmemory: bool,
    nthreads: usize,
) -> Result<(), BBIReadError> {
    let runtime = runtime::Builder::new_multi_thread()
        .worker_threads(nthreads)
        .build()
        .unwrap();

    let mut remaining_chroms = bigwig.chroms().to_vec();
    remaining_chroms.reverse();

    let mut chrom_files: VecDeque<_> = VecDeque::new();

    async fn file_future<R: SeekableRead + 'static>(
        mut bigwig: BigWigRead<R>,
        chrom: ChromInfo,
        mut writer: io::BufWriter<TempFileBufferWriter<File>>,
    ) -> Result<(), BBIReadError> {
        let mut buf: String = String::with_capacity(50); // Estimate
        for raw_val in bigwig.get_interval(&chrom.name, 0, chrom.length)? {
            let val = raw_val?;
            // Using ryu for f32 to string conversion has a ~15% speedup
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
            buf.clear();
        }
        Ok(())
    }

    loop {
        while chrom_files.len() < nthreads {
            let Some(chrom) = remaining_chroms.pop() else {
                break;
            };

            let bigbed = bigwig.reopen()?;
            let (buf, file): (TempFileBuffer<File>, TempFileBufferWriter<File>) =
                TempFileBuffer::new(inmemory);
            let writer = io::BufWriter::new(file);
            let handle = runtime
                .spawn(file_future(bigbed, chrom, writer))
                .map(|f| f.unwrap());
            chrom_files.push_back((handle, buf));
        }

        let Some((f, mut buf)) = chrom_files.pop_front() else {
            break;
        };

        buf.switch(out_file);
        runtime.block_on(f).unwrap();
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
