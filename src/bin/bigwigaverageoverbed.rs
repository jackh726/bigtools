use std::collections::VecDeque;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Seek, SeekFrom, Write};

use bigtools::bed::bedparser::{parse_bed, BedParser};
use bigtools::bed::indexer::index_chroms;
use bigtools::utils::chromvalues::ChromValues;
use bigtools::utils::reopen::{Reopen, SeekableRead};
use bigtools::utils::streaming_linereader::StreamingLineReader;
use clap::{Arg, Command};

use bigtools::bbi::BigWigRead;
use bigtools::utils::misc::{stats_for_bed_item, Name};
use crossbeam_channel::TryRecvError;

fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
    let matches = Command::new("BigWigAverageOverBed")
        .arg(Arg::new("bigwig")
                .help("The input bigwig file")
                .index(1)
                .required(true)
            )
        .arg(Arg::new("bedin")
                .help("The input bed file")
                .index(2)
                .required(true)
            )
        .arg(Arg::new("output")
                .help("The output bed file")
                .index(3)
                .required(true)
            )
        .arg(Arg::new("namecol")
                .short('n')
                .help("Supports three types of options: `interval`, `none`, or a column number (one indexed). If `interval`, the name column in the output will be the interval in the form of `chrom:start-end`. If `none`, then all columns will be included in the output file. Otherwise, the one-indexed column will be used as the name. By default, column 4 is used as a name column.")
                .default_value("4")
            )
        .arg(Arg::new("nthreads")
            .short('t')
            .help("Number of threads to use. Defaults to 1.")
            .default_value("1")
            )
        .get_matches();

    let bigwigpath = matches.get_one::<String>("bigwig").unwrap();
    let bedinpath = matches.get_one::<String>("bedin").unwrap();
    let bedoutpath = matches.get_one::<String>("output").unwrap();

    let mut inbigwig = BigWigRead::open_file(bigwigpath)?;
    let outbed = File::create(bedoutpath)?;
    let mut bedoutwriter = BufWriter::new(outbed);

    let bedin = BufReader::new(File::open(bedinpath)?);
    let mut bedstream = StreamingLineReader::new(bedin);

    let name = match matches.get_one::<String>("namecol").map(String::as_ref) {
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

    let nthreads: usize = *matches.get_one::<usize>("nthreads").unwrap();
    let parallel = nthreads > 1;

    if parallel {
        fn process_chrom<R: Reopen + SeekableRead>(
            start: u64,
            chrom: String,
            bedinpath: String,
            name: Name,
            inbigwig: &mut BigWigRead<R>,
        ) -> Result<File, Box<dyn Error + Send + Sync>> {
            let mut tmp = tempfile::tempfile()?;

            let mut chrom_bed_file = File::open(bedinpath)?;
            chrom_bed_file.seek(SeekFrom::Start(start))?;
            let mut bed_parser = BedParser::from_bed_file(chrom_bed_file);
            let mut data = match bed_parser.next_chrom() {
                Some(Ok((next_chrom, data))) => {
                    if next_chrom != chrom {
                        return Err(io::Error::new(
                            io::ErrorKind::Other,
                            format!("Invalid indexing."),
                        )
                        .into());
                    }
                    data
                }
                Some(Err(e)) => {
                    return Err(e.into());
                }
                None => {
                    return Err(
                        io::Error::new(io::ErrorKind::Other, format!("Invalid indexing.")).into(),
                    );
                }
            };

            loop {
                let entry = match data.next() {
                    None => break,
                    Some(Err(e)) => {
                        return Err(e.into());
                    }
                    Some(Ok(entry)) => entry,
                };

                let entry = match stats_for_bed_item(name, &chrom, entry, inbigwig) {
                    Ok(stats) => stats,
                    Err(e) => {
                        return Err(e.into());
                    }
                };

                let stats = format!(
                    "{}\t{}\t{:.3}\t{:.3}\t{:.3}",
                    entry.size, entry.bases, entry.sum, entry.mean0, entry.mean
                );
                writeln!(&mut tmp, "{}\t{}", entry.name, stats)?
            }

            Ok(tmp)
        }

        let bed = File::open(bedinpath)?;
        let chrom_indices: Vec<(u64, String)> = index_chroms(bed)?;

        let mut chrom_data = VecDeque::with_capacity(chrom_indices.len());
        let (chrom_data_sender, chrom_data_receiver) = crossbeam_channel::unbounded();
        for (start, chrom) in chrom_indices {
            let (result_sender, result_receiver) = crossbeam_channel::bounded(1);
            chrom_data.push_back(result_receiver);
            chrom_data_sender
                .send((start, chrom, bedinpath.to_string(), result_sender))
                .unwrap();
        }
        drop(chrom_data_sender);
        let mut threads = Vec::with_capacity(nthreads - 1);
        for _ in 0..(nthreads - 1) {
            let inbigwig_ = inbigwig.reopen()?;
            let chrom_data_receiver_ = chrom_data_receiver.clone();
            let do_process_chrom = move || {
                let mut inbigwig = inbigwig_;
                let chrom_data_receiver = chrom_data_receiver_;
                loop {
                    let next_chrom = chrom_data_receiver.recv();
                    let (start, chrom, bedinpath, result_sender) = match next_chrom {
                        Ok(n) => n,
                        Err(_) => break,
                    };

                    let result = process_chrom(start, chrom, bedinpath, name, &mut inbigwig);
                    result_sender.send(result).unwrap();
                }
            };
            let join_handle = std::thread::spawn(do_process_chrom);
            threads.push(join_handle);
        }
        while let Some(result_receiver) = chrom_data.pop_front() {
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
                            let next_chrom = chrom_data_receiver.recv();
                            let (start, chrom, bedinpath, result_sender) = match next_chrom {
                                Ok(n) => n,
                                Err(_) => {
                                    wait = true;
                                    continue;
                                }
                            };

                            let result =
                                process_chrom(start, chrom, bedinpath, name, &mut inbigwig);
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

            let entry = match stats_for_bed_item(name, chrom, entry, &mut inbigwig) {
                Ok(stats) => stats,
                Err(e) => return Err(e.into()),
            };

            let stats = format!(
                "{}\t{}\t{:.3}\t{:.3}\t{:.3}",
                entry.size, entry.bases, entry.sum, entry.mean0, entry.mean
            );
            writeln!(&mut bedoutwriter, "{}\t{}", entry.name, stats)?
        }
    }

    Ok(())
}
