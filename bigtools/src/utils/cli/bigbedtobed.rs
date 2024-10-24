use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::Path;

use clap::Parser;
use futures::{SinkExt, StreamExt};
use tokio::runtime;
use ufmt::uwrite;

use crate::utils::reopen::{Reopen, SeekableRead};
use crate::utils::streaming_linereader::StreamingLineReader;
use crate::utils::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};
use crate::{BBIReadError, BigBedRead, ChromInfo};

#[derive(Clone, Debug, PartialEq, Parser)]
#[command(
    name = "bigbedtobed",
    about = "Converts an input bigBed to a bed.",
    long_about = "Converts an input bigBed to a bed. Can be multi-threaded for substantial speedups. Note for roughly each core, one temporary file will be opened."
)]
pub struct BigBedToBedArgs {
    /// the bigbed to get convert to bed
    pub big_bed: String,

    /// the path of the bed to output to
    pub bed: String,

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

    /// Set the number of threads to use. This tool will nearly always benefit from more cores (<= # chroms). Note: for parts of the runtime, the actual usage may be nthreads+1
    #[arg(short = 't', long)]
    #[arg(default_value_t = 6)]
    pub nthreads: usize,

    /// Do not create temporary files for intermediate data. (Only applicable when using multiple threads.)
    /// By default, approximately one temporary file will be opened for each core.
    #[arg(long)]
    #[arg(default_value_t = false)]
    pub inmemory: bool,

    /// If set, outputs the values for a given zoom level. The format will be
    /// <chrom> <zoom_start>    <zoom_end>  <total_items> <bases_covered>   <min_val>   <max_val>   <sum>   <sum_squares>
    #[arg(long)]
    pub zoom: Option<u32>,
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
            // Right now, we don't offload decompression to separate threads,
            // so specifying `chrom` effectively means single-threaded
            if nthreads == 1 || args.chrom.is_some() || args.zoom.is_some() {
                write_bed_singlethreaded(bigbed, bed, args.chrom, args.start, args.end, args.zoom)?;
            } else {
                write_bed(bigbed, bed, args.inmemory, nthreads)?;
            }
        }
    }

    Ok(())
}

pub fn write_bed_singlethreaded<R: Reopen + SeekableRead>(
    mut bigbed: BigBedRead<R>,
    out_file: File,
    chrom: Option<String>,
    start: Option<u32>,
    end: Option<u32>,
    zoom: Option<u32>,
) -> Result<(), Box<dyn Error>> {
    let start = chrom.as_ref().and_then(|_| start);
    let end = chrom.as_ref().and_then(|_| end);

    let chroms: Vec<ChromInfo> = if let Some(arg_chrom) = chrom {
        let chrom = bigbed.chroms().iter().find(|c| c.name == arg_chrom);
        let Some(chrom) = chrom else {
            eprintln!("{arg_chrom} not found in file.");
            return Ok(());
        };
        vec![chrom.clone()]
    } else {
        bigbed.chroms().to_vec()
    };
    let mut writer = io::BufWriter::with_capacity(32 * 1000, out_file);
    let mut buf: String = String::with_capacity(50); // Estimate
    if let Some(zoom) = zoom {
        for chrom in chroms {
            let start = start.unwrap_or(0);
            let end = end.unwrap_or(chrom.length);
            for raw_val in bigbed
                .get_zoom_interval(&chrom.name, start, end, zoom)
                .unwrap()
            {
                let val = raw_val?;
                uwrite!(
                    &mut buf,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    chrom.name,
                    val.start,
                    val.end,
                    val.summary.total_items,
                    val.summary.bases_covered,
                    val.summary.min_val.to_string(),
                    val.summary.max_val.to_string(),
                    val.summary.sum.to_string(),
                    val.summary.sum_squares.to_string()
                )?;
                writer.write(buf.as_bytes())?;
                buf.clear();
            }
        }
    } else {
        for chrom in chroms {
            let start = start.unwrap_or(0);
            let end = end.unwrap_or(chrom.length);
            for raw_val in bigbed.get_interval(&chrom.name, start, end)? {
                let val = raw_val?;
                if !val.rest.is_empty() {
                    uwrite!(
                        &mut buf,
                        "{}\t{}\t{}\t{}\n",
                        chrom.name,
                        val.start,
                        val.end,
                        val.rest
                    )
                    .unwrap();
                } else {
                    uwrite!(&mut buf, "{}\t{}\t{}\n", chrom.name, val.start, val.end).unwrap();
                };
                writer.write(buf.as_bytes())?;
                buf.clear();
            }
        }
    }
    Ok(())
}

pub fn write_bed<R: Reopen + SeekableRead + Send + 'static>(
    bigbed: BigBedRead<R>,
    mut out_file: File,
    inmemory: bool,
    nthreads: usize,
) -> Result<(), BBIReadError> {
    let runtime = runtime::Builder::new_multi_thread()
        .worker_threads(nthreads)
        .build()
        .unwrap();

    let mut remaining_chroms = bigbed.chroms().to_vec();
    remaining_chroms.reverse();

    async fn file_future<R: SeekableRead + 'static>(
        mut bigbed: BigBedRead<R>,
        chrom: ChromInfo,
        mut writer: io::BufWriter<TempFileBufferWriter<File>>,
    ) -> Result<(), BBIReadError> {
        let mut buf: String = String::with_capacity(50); // Estimate
        for raw_val in bigbed.get_interval(&chrom.name, 0, chrom.length)? {
            let val = raw_val?;
            if !val.rest.is_empty() {
                uwrite!(
                    &mut buf,
                    "{}\t{}\t{}\t{}\n",
                    chrom.name,
                    val.start,
                    val.end,
                    val.rest
                )
                .unwrap();
            } else {
                uwrite!(&mut buf, "{}\t{}\t{}\n", chrom.name, val.start, val.end).unwrap();
            };
            writer.write(buf.as_bytes())?;
            buf.clear();
        }
        Ok(())
    }

    let (mut handle_snd, mut handle_rcv) = futures::channel::mpsc::channel(nthreads);
    let (mut buf_snd, mut buf_rcv) = futures::channel::mpsc::unbounded();
    runtime.spawn(async move {
        loop {
            let Some(chrom) = remaining_chroms.pop() else {
                return Ok::<_, BBIReadError>(());
            };

            let bigbed = bigbed.reopen()?;
            let (buf, file): (TempFileBuffer<File>, TempFileBufferWriter<File>) =
                TempFileBuffer::new(inmemory);
            let writer = io::BufWriter::new(file);
            let handle = tokio::task::spawn(file_future(bigbed, chrom, writer));

            handle_snd.send(handle).await.unwrap();
            buf_snd.send(buf).await.unwrap();
        }
    });

    let data_handle = runtime.spawn(async move {
        loop {
            let next = handle_rcv.next().await;
            let Some(handle) = next else {
                return Ok::<_, BBIReadError>(());
            };
            handle.await.unwrap()?;
        }
    });
    runtime.block_on(async move {
        loop {
            let next = buf_rcv.next().await;
            let Some(mut buf) = next else {
                data_handle.await.unwrap()?;
                return Ok::<_, BBIReadError>(());
            };

            buf.switch(out_file);
            while !buf.is_real_file_ready() {
                tokio::task::yield_now().await;
            }
            out_file = buf.await_real_file();
        }
    })?;

    Ok(())
}

pub fn write_bed_from_bed<R: Reopen + SeekableRead + Send + 'static>(
    mut bigbed: BigBedRead<R>,
    out_file: File,
    bed: File,
) -> Result<(), BBIReadError> {
    let mut bedstream = StreamingLineReader::new(BufReader::new(bed));
    let mut writer = io::BufWriter::new(out_file);

    let mut buf = String::with_capacity(50); // Estimate
    while let Some(line) = bedstream.read() {
        let line = line?;
        let mut split = line.trim().splitn(5, '\t');
        let chrom = split.next().expect("Missing chrom");
        let start = split.next().expect("Missing start").parse::<u32>().unwrap();
        let end = split.next().expect("Missing end").parse::<u32>().unwrap();
        for raw_val in bigbed.get_interval(chrom, start, end)? {
            let mut val = raw_val?;
            val.start = val.start.max(start);
            val.end = val.end.min(end);
            if !val.rest.is_empty() {
                uwrite!(
                    &mut buf,
                    "{}\t{}\t{}\t{}\n",
                    chrom,
                    val.start,
                    val.end,
                    val.rest
                )
                .unwrap();
            } else {
                uwrite!(&mut buf, "{}\t{}\t{}\n", chrom, val.start, val.end).unwrap();
            };
            writer.write(buf.as_bytes())?;
            buf.clear();
        }
    }

    Ok(())
}
