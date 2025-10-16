use std::collections::VecDeque;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Seek, SeekFrom, Write};

use clap::Parser;
use crossbeam_channel::TryRecvError;

use crate::bed::bedparser::{parse_bed, BedFileStream, StreamingBedValues};
use crate::utils::file_view::FileView;
use crate::utils::misc::{name_for_bed_item, stats_for_bed_item, BigWigAverageOverBedEntry, Name};
use crate::utils::reopen::{Reopen, ReopenableFile};
use crate::utils::split_file_into_chunks_by_size;
use crate::utils::streaming_linereader::StreamingLineReader;
use crate::{BBIFileRead, BBIReadError, BigWigRead};

#[derive(Clone, Debug, PartialEq, Parser)]
#[command(
    name = "bigwigaverageoverbed",
    about = "Gets statistics of a bigWig over a bed.",
    long_about = None,
)]
pub struct BigWigAverageOverBedArgs {
    /// The input bigwig file
    pub bigwig: String,

    /// The input bed file.
    pub bedin: String,

    /// The output bed file.
    /// Specifying `-` will write to `stdout`.
    pub output: String,

    /// Supports three types of options: `interval`, `none`, or a column number (one indexed).
    /// If `interval`, the name column in the output will be the interval in the form of `chrom:start-end`.
    /// If `none`, then all columns will be included in the output file.
    /// Otherwise, the one-indexed column will be used as the name. By default, column 4 is used as a name column.
    #[arg(short = 'n', long)]
    pub namecol: Option<String>,

    /// If set, restrict output to given chromosome
    #[arg(long)]
    pub chrom: Option<String>,

    /// If set, restrict output to regions greater than or equal to it
    #[arg(long)]
    pub start: Option<u32>,

    /// If set, restrict output to regions less than it
    #[arg(long)]
    pub end: Option<u32>,

    /// If set, include two additional columns containing the min and max observed in the area
    #[arg(long)]
    pub min_max: bool,

    /// Set the number of threads to use. This tool will nearly always benefit from more cores.
    /// Note: for parts of the runtime, the actual usage may be nthreads+1
    #[arg(short = 't', long)]
    #[arg(default_value_t = 6)]
    pub nthreads: usize,
}

fn bigwigaverageoverbed_impl<O: Write>(
    args: BigWigAverageOverBedArgs,
    bedoutwriter: O,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let bigwigpath = args.bigwig;
    let bedinpath = args.bedin;

    let add_min_max = args.min_max;

    let reopen = ReopenableFile {
        file: File::open(&bigwigpath)?,
        path: bigwigpath.into(),
    };
    let mut inbigwig = BigWigRead::open(reopen)?.cached();

    let mut bedoutwriter = BufWriter::new(bedoutwriter);

    let name = match args.namecol.as_deref() {
        Some("interval") => Name::Interval,
        Some("none") => Name::None,
        Some(col) => {
            let col = col.parse::<usize>();
            match col {
                Ok(col) if col > 0 => Name::Column(col - 1),
                Ok(_) => return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid name column option. Column values are one-indexed, so should not be zero.",
                ).into()),
                Err(_) => return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid name column option. Allowed options are `interval`, `none`, or an one-indexed integer value for a given column.",
                ).into()),
            }
        }
        None => Name::Column(3),
    };

    let nthreads: usize = args.nthreads;
    let parallel = nthreads > 1;

    if parallel {
        let chunks = split_file_into_chunks_by_size(File::open(&bedinpath)?, nthreads as u64)?;

        fn process_chunk<R: Reopen + BBIFileRead>(
            start: u64,
            end: u64,
            bedinpath: String,
            name: Name,
            add_min_max: bool,
            inbigwig: &mut BigWigRead<R>,
        ) -> Result<File, Box<dyn Error + Send + Sync>> {
            let mut tmp = tempfile::tempfile()?;

            let chrom_bed_file = File::open(bedinpath)?;
            let chrom_bed_file = FileView::new(chrom_bed_file, start, end)?;
            let mut bed_stream = BedFileStream {
                bed: StreamingLineReader::new(BufReader::new(chrom_bed_file)),
                parse: parse_bed,
            };

            loop {
                let (chrom, entry) = match bed_stream.next() {
                    None => break,
                    Some(Err(e)) => {
                        return Err(e.into());
                    }
                    Some(Ok(entry)) => entry,
                };

                let name = match name_for_bed_item(name, chrom, &entry) {
                    Ok(name) => name,
                    Err(e) => {
                        return Err(e.into());
                    }
                };

                let size = entry.end - entry.start;

                let entry = match stats_for_bed_item(chrom, entry, inbigwig) {
                    Ok(stats) => stats,
                    Err(BBIReadError::InvalidChromosome(..)) => BigWigAverageOverBedEntry {
                        bases: 0,
                        max: 0.0,
                        min: 0.0,
                        mean: 0.0,
                        mean0: 0.0,
                        size,
                        sum: 0.0,
                    },
                    Err(e) => {
                        return Err(e.into());
                    }
                };

                let stats = match add_min_max {
                    true => format!(
                        "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                        entry.size,
                        entry.bases,
                        entry.sum,
                        entry.mean0,
                        entry.mean,
                        entry.min,
                        entry.max
                    ),
                    false => format!(
                        "{}\t{}\t{:.3}\t{:.3}\t{:.3}",
                        entry.size, entry.bases, entry.sum, entry.mean0, entry.mean
                    ),
                };
                writeln!(&mut tmp, "{}\t{}", name, stats)?
            }

            Ok(tmp)
        }

        let mut chunk_data = VecDeque::with_capacity(chunks.len());
        let (chunk_data_sender, chunk_data_receiver) = crossbeam_channel::unbounded();
        for (start, end) in chunks {
            let (result_sender, result_receiver) = crossbeam_channel::bounded(1);
            chunk_data.push_back(result_receiver);
            chunk_data_sender
                .send((start, end, bedinpath.to_string(), result_sender))
                .unwrap();
        }
        drop(chunk_data_sender);
        let mut threads = Vec::with_capacity(nthreads);
        for _ in 0..(nthreads) {
            let inbigwig_ = inbigwig.reopen()?;
            let chunk_data_receiver_ = chunk_data_receiver.clone();
            let do_process_chrom = move || {
                let mut inbigwig = inbigwig_;
                let chunk_data_receiver = chunk_data_receiver_;
                loop {
                    let next_chunk = chunk_data_receiver.recv();
                    let (start, end, bedinpath, result_sender) = match next_chunk {
                        Ok(n) => n,
                        Err(_) => break,
                    };

                    let result =
                        process_chunk(start, end, bedinpath, name, add_min_max, &mut inbigwig);
                    result_sender.send(result).unwrap();
                }
            };
            let join_handle = std::thread::spawn(do_process_chrom);
            threads.push(join_handle);
        }
        while let Some(result_receiver) = chunk_data.pop_front() {
            let mut wait = false;
            loop {
                if !wait {
                    let result = result_receiver.try_recv();
                    match result {
                        Ok(result) => {
                            let mut tmp = result?;
                            tmp.seek(SeekFrom::Start(0))?;
                            io::copy(&mut tmp, &mut bedoutwriter)?;
                            break;
                        }
                        Err(e @ TryRecvError::Disconnected) => return Err(e.into()),
                        Err(TryRecvError::Empty) => {
                            let next_chunk = chunk_data_receiver.recv();
                            let (start, chrom, bedinpath, result_sender) = match next_chunk {
                                Ok(n) => n,
                                Err(_) => {
                                    wait = true;
                                    continue;
                                }
                            };

                            let result = process_chunk(
                                start,
                                chrom,
                                bedinpath,
                                name,
                                add_min_max,
                                &mut inbigwig,
                            );
                            result_sender.send(result).unwrap();
                        }
                    }
                } else {
                    let result = result_receiver.recv();
                    match result {
                        Ok(result) => {
                            let mut tmp = result?;
                            tmp.seek(SeekFrom::Start(0))?;
                            io::copy(&mut tmp, &mut bedoutwriter)?;
                            break;
                        }
                        Err(e) => return Err(e.into()),
                    }
                }
            }
        }
    } else {
        let bed = File::open(&bedinpath)?;
        let bedin: BufReader<File> = BufReader::new(bed);
        let mut bedstream = StreamingLineReader::new(bedin);

        loop {
            let line = match bedstream.read() {
                None => break,
                Some(Err(e)) => {
                    return Err(e.into());
                }
                Some(Ok(line)) => line,
            };

            let (chrom, entry) = parse_bed(line).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid bed: A minimum of 3 columns must be specified (chrom, start, end).",
                )
            })??;

            let name = match name_for_bed_item(name, chrom, &entry) {
                Ok(name) => name,
                Err(e) => {
                    return Err(e.into());
                }
            };

            let entry = match stats_for_bed_item(chrom, entry, &mut inbigwig) {
                Ok(stats) => stats,
                Err(e) => return Err(e.into()),
            };

            let stats = match add_min_max {
                true => format!(
                    "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                    entry.size,
                    entry.bases,
                    entry.sum,
                    entry.mean0,
                    entry.mean,
                    entry.min,
                    entry.max
                ),
                false => format!(
                    "{}\t{}\t{:.3}\t{:.3}\t{:.3}",
                    entry.size, entry.bases, entry.sum, entry.mean0, entry.mean
                ),
            };
            writeln!(&mut bedoutwriter, "{}\t{}", name, stats)?
        }
    }

    Ok(())
}

pub fn bigwigaverageoverbed(
    args: BigWigAverageOverBedArgs,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let bedoutpath = args.output.clone();
    match bedoutpath.as_str() {
        "-" => {
            let stdout = io::stdout().lock();
            bigwigaverageoverbed_impl(args, stdout)
        }
        _ => {
            let outbed = File::create(bedoutpath)?;
            bigwigaverageoverbed_impl(args, outbed)
        }
    }
}
