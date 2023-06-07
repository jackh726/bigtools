use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use clap::{App, Arg};
use thiserror::Error;

use bigtools::bbi::Value;
use bigtools::bbi::{BBIRead, BigWigRead, BigWigWrite};
use bigtools::bbiread::BBIReadError;
use bigtools::utils::chromvalues::ChromValues;
use bigtools::utils::filebufferedchannel;
use bigtools::utils::merge::merge_sections_many;
use bigtools::utils::seekableread::ReopenableFile;
use bigtools::{ChromData, ChromDataState, ChromProcessingFnOutput};

pub struct MergingValues {
    // We Box<dyn Iterator> because other this would be a mess to try to type
    iter: std::iter::Peekable<Box<dyn Iterator<Item = Result<Value, MergingValuesError>> + Send>>,
}

impl MergingValues {
    pub fn new<I: 'static>(iters: Vec<I>) -> Self
    where
        I: Iterator<Item = Result<Value, MergingValuesError>> + Send,
    {
        let iter: Box<dyn Iterator<Item = Result<Value, MergingValuesError>> + Send> = Box::new(
            merge_sections_many(iters)
                .filter(|x| x.as_ref().map(|v| v.value != 0.0).unwrap_or(true)),
        );
        MergingValues {
            iter: iter.peekable(),
        }
    }
}

#[derive(Error, Debug)]
pub enum MergingValuesError {
    #[error("{}", .0)]
    BBIReadError(#[from] BBIReadError),
    #[error("{}", .0)]
    MismatchedChroms(String),
    #[error("{}", .0)]
    Other(String),
    #[error("{}", .0)]
    IoError(#[from] io::Error),
}

impl ChromValues for MergingValues {
    type Value = Value;
    type Error = MergingValuesError;

    fn next(&mut self) -> Option<Result<Value, MergingValuesError>> {
        match self.iter.next() {
            Some(Ok(v)) => Some(Ok(v)),
            Some(Err(e)) => Some(Err(e.into())),
            None => None,
        }
    }

    fn peek(&mut self) -> Option<Result<&Value, &MergingValuesError>> {
        match self.iter.peek() {
            Some(Ok(v)) => Some(Ok(v)),
            Some(Err(err)) => Some(Err(err)),
            None => None,
        }
    }
}

pub fn get_merged_vals(
    bigwigs: Vec<BigWigRead<ReopenableFile, File>>,
    max_zooms: usize,
) -> Result<
    (
        impl Iterator<Item = Result<(String, u32, MergingValues), MergingValuesError>>,
        HashMap<String, u32>,
    ),
    MergingValuesError,
> {
    let (chrom_sizes, chrom_map) = {
        // NOTE: We don't need to worry about max fds here because chroms are cached.

        // Get sizes for each and check that all files (that have the chrom) agree
        // Check that all chrom sizes match for all files
        let mut chrom_sizes = BTreeMap::new();
        let mut chrom_map = HashMap::new();
        for chrom in bigwigs.iter().flat_map(BBIRead::get_chroms).map(|c| c.name) {
            if chrom_sizes.get(&chrom).is_some() {
                continue;
            }
            let (sizes, bws): (Vec<u32>, Vec<BigWigRead<ReopenableFile, File>>) = bigwigs
                .iter()
                .map(|w| {
                    let chroms = w.get_chroms();
                    let res = chroms.iter().find(|v| v.name == chrom);
                    res.map(|s| (s.length, w.clone()))
                })
                .flatten()
                .unzip();
            let size = sizes[0];
            if !sizes.iter().all(|s| *s == size) {
                eprintln!("Chrom '{:?}' had different sizes in the bigwig files. (Are you using the same assembly?)", chrom);
                return Err(MergingValuesError::MismatchedChroms(
                    "Invalid input (nonmatching chroms)".to_owned(),
                ));
            }

            chrom_sizes.insert(chrom.clone(), (size, bws));
            chrom_map.insert(chrom.clone(), size);
        }

        (chrom_sizes, chrom_map)
    };

    const MAX_FDS: usize = 1000;
    const PARALLEL_CHROMS: usize = 1;
    // This might be a *bit* conservative, but is really mostly an estimate
    let max_bw_fds: usize = MAX_FDS
        - 1 /* output bigWig (data) */
        - 1 /* index */
        - (1 /* data sections */ + 1  /* index sections */ + max_zooms /* zoom data sections */ + max_zooms /* zoom index sections */) * PARALLEL_CHROMS;

    let iter = chrom_sizes.into_iter().map(move |(chrom, (size, bws))| {
        if bws.len() > max_bw_fds {
            eprintln!("Number of bigWigs to merge would exceed the maximum number of file descriptors. Splitting into chunks.");

            let mut merges: Vec<Box<dyn Iterator<Item = Result<Value, MergingValuesError>> + Send>> = bws
                .into_iter()
                .map(|b| {
                    let iter = b.get_interval_move(&chrom, 1, size).map(|i| i.map(|r| r.map_err(|e| MergingValuesError::BBIReadError(e))))?;
                    Ok(Box::new(iter) as Box<_>)
                })
                .collect::<Result<Vec<_>, BBIReadError>>()?;

            while merges.len() > max_bw_fds {
                merges = {
                    let len = merges.len();
                    let mut vals = merges.into_iter().peekable();
                    let mut merges: Vec<Box<dyn Iterator<Item = Result<Value, MergingValuesError>> + Send>> = Vec::with_capacity(len/max_bw_fds+1);

                    while vals.peek().is_some() {
                        let chunk = vals.by_ref().take(max_bw_fds).collect::<Vec<_>>();
                        let mut mergingvalues = MergingValues::new(chunk);
                        let (mut sender, receiver) = filebufferedchannel::lazy_channel::<Value>(3200)?;
                        while let Some(val) = mergingvalues.next() {
                            let val = val?;
                            sender.send(val).unwrap();
                        }

                        merges.push(Box::new(receiver.into_iter().map(Result::Ok)));
                    }
                    merges
                };
            }

            let mergingvalues = MergingValues::new(merges);
            Ok((chrom, size, mergingvalues))
        } else {
            let iters: Vec<_> = bws
                .into_iter()
                .map(|b| b.get_interval_move(&chrom, 1, size).map(|i| i.map(|r| r.map_err(|e| MergingValuesError::BBIReadError(e)))))
                .collect::<Result<Vec<_>, _>>()?;
            let mergingvalues = MergingValues::new(iters);

            Ok((chrom, size, mergingvalues))
        }
    });

    Ok((iter, chrom_map))
}

struct ChromGroupReadImpl {
    iter: Box<dyn Iterator<Item = Result<(String, u32, MergingValues), MergingValuesError>> + Send>,
}

impl<E: From<io::Error>> ChromData<E> for ChromGroupReadImpl {
    type Output = MergingValues;

    fn advance<
        F: FnMut(String, Self::Output) -> Result<ChromProcessingFnOutput<Self::Output>, E>,
    >(
        &mut self,
        do_read: &mut F,
    ) -> Result<ChromDataState<Self::Output>, E> {
        let next = self.iter.next();
        Ok(match next {
            Some(Err(err)) => ChromDataState::Error(err.into()),
            Some(Ok((chrom, _, mergingvalues))) => {
                let read = do_read(chrom, mergingvalues)?;

                ChromDataState::NewChrom(read)
            }
            None => ChromDataState::Finished,
        })
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("BigWigMerge")
        .arg(Arg::new("output")
                .help("the path of the merged output bigwig (if .bw or .bigWig) or bedGraph (if .bedGraph)")
                .index(1)
                .required(true)
            )
        .arg(Arg::new("bigwig")
                .short('b')
                .help("the path of an input bigwig to merge")
                .multiple_occurrences(true)
                .takes_value(true)
            )
        .arg(Arg::new("list")
                .short('l')
                .help("a line-delimited list of bigwigs")
                .multiple_occurrences(true)
                .takes_value(true)
            )
        .arg(Arg::new("nthreads")
                .short('t')
                .help("Set the number of threads to use")
                .takes_value(true)
                .default_value("6"))
        .get_matches();

    let output = matches.value_of("output").unwrap().to_owned();
    let mut bigwigs: Vec<BigWigRead<ReopenableFile, File>> = vec![];

    if let Some(bws) = matches.values_of("bigwig") {
        for name in bws {
            match BigWigRead::from_file_and_attach(name).map(|mut bw| {
                bw.close();
                bw
            }) {
                Ok(bw) => bigwigs.push(bw),
                Err(e) => {
                    eprintln!("Error when opening bigwig ({}): {:?}", name, e);
                    return Ok(());
                }
            }
        }
    }
    if let Some(lists) = matches.values_of("list") {
        for list in lists {
            let list_file = match File::open(list) {
                Ok(f) => f,
                Err(e) => {
                    eprintln!("Couldn't open file: {:?}", e);
                    return Ok(());
                }
            };
            let lines = BufReader::new(list_file).lines();
            for line in lines {
                let name = line?;
                match BigWigRead::from_file_and_attach(&name).map(|mut bw| {
                    bw.close();
                    bw
                }) {
                    Ok(bw) => bigwigs.push(bw),
                    Err(e) => {
                        eprintln!("Error when opening bigwig ({}): {:?}", name, e);
                        return Ok(());
                    }
                }
            }
        }
    }

    let nthreads = {
        let nthreads = matches.value_of("nthreads").unwrap();
        let parsed = nthreads.parse();
        if parsed.is_err() {
            eprintln!("Invalid argument for `nthreads`: must be a positive number");
            return Ok(());
        }
        parsed.unwrap()
    };

    let (iter, chrom_map) = get_merged_vals(bigwigs, 10)?;

    match output {
        output if output.ends_with(".bw") || output.ends_with(".bigWig") => {
            let outb = BigWigWrite::create_file(output);
            let pool = futures::executor::ThreadPoolBuilder::new()
                .pool_size(nthreads)
                .create()
                .expect("Unable to create thread pool.");
            let all_values = ChromGroupReadImpl {
                iter: Box::new(iter),
            };
            outb.write(chrom_map, all_values, pool)?;
        }
        output if output.ends_with(".bedGraph") => {
            // TODO: convert to multi-threaded
            use std::io::Write;

            let bedgraph = File::create(output)?;
            let mut writer = io::BufWriter::new(bedgraph);

            for v in iter {
                let (chrom, _, mut values) = v?;
                while let Some(val) = values.next() {
                    let val = val?;
                    writer.write_fmt(format_args!(
                        "{}\t{}\t{}\t{}\n",
                        chrom, val.start, val.end, val.value
                    ))?;
                }
            }
        }
        _ => {
            eprintln!("Invalid output file. Must end with .bw or .bigWig for bigwig or .bedGraph for bedGraph");
            return Ok(());
        }
    }

    //TODO: fails with too many open files
    Ok(())
}
