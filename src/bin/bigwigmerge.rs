use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use bigtools::{ChromData, ChromDataState};
use clap::{App, Arg};

use bigtools::bigwig::BBIWriteOptions;
use bigtools::bigwig::Value;
use bigtools::bigwig::{BBIRead, BigWigRead, BigWigWrite, WriteGroupsError};
use bigtools::utils::chromvalues::ChromValues;
use bigtools::utils::filebufferedchannel;
use bigtools::utils::idmap::IdMap;
use bigtools::utils::merge::merge_sections_many;
use bigtools::utils::seekableread::ReopenableFile;

pub struct MergingValues {
    // We Box<dyn Iterator> because other this would be a mess to try to type
    iter: std::iter::Peekable<Box<dyn Iterator<Item = io::Result<Value>> + Send>>,
}

impl MergingValues {
    pub fn new<I: 'static>(iters: Vec<I>) -> Self
    where
        I: Iterator<Item = io::Result<Value>> + Send,
    {
        let iter: Box<dyn Iterator<Item = io::Result<Value>> + Send> = Box::new(
            merge_sections_many(iters)
                .filter(|x| x.as_ref().map(|v| v.value != 0.0).unwrap_or(true)),
        );
        MergingValues {
            iter: iter.peekable(),
        }
    }
}

impl ChromValues<Value> for MergingValues {
    fn next(&mut self) -> Option<io::Result<Value>> {
        self.iter.next()
    }

    fn peek(&mut self) -> Option<&Value> {
        match self.iter.peek() {
            Some(Ok(v)) => Some(v),
            Some(Err(_)) => None,
            None => None,
        }
    }
}

pub fn get_merged_vals(
    bigwigs: Vec<BigWigRead<ReopenableFile, File>>,
    max_zooms: usize,
) -> io::Result<(
    impl Iterator<Item = io::Result<(String, u32, MergingValues)>>,
    HashMap<String, u32>,
)> {
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
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "Invalid input (nonmatching chroms)",
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

            let mut merges: Vec<Box<dyn Iterator<Item = io::Result<Value>> + Send>> = bws
                .into_iter()
                .map(|b| {
                    let iter = b.get_interval_move(&chrom, 1, size)?;
                    Ok(Box::new(iter) as Box<_>)
                })
                .collect::<io::Result<Vec<_>>>()?;

            while merges.len() > max_bw_fds {
                merges = {
                    let len = merges.len();
                    let mut vals = merges.into_iter().peekable();
                    let mut merges: Vec<Box<dyn Iterator<Item = io::Result<Value>> + Send>> = Vec::with_capacity(len/max_bw_fds+1);

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
                .map(|b| b.get_interval_move(&chrom, 1, size))
                .collect::<io::Result<Vec<_>>>()?;
            let mergingvalues = MergingValues::new(iters);

            Ok((chrom, size, mergingvalues))
        }
    });

    Ok((iter, chrom_map))
}

struct ChromGroupReadImpl {
    pool: futures::executor::ThreadPool,
    options: BBIWriteOptions,
    iter: Box<dyn Iterator<Item = io::Result<(String, u32, MergingValues)>> + Send>,
    chrom_ids: Option<IdMap>,
}

impl ChromData for ChromGroupReadImpl {
    fn advance(mut self) -> ChromDataState<Self> {
        let next = self.iter.next();
        match next {
            Some(Err(err)) => ChromDataState::Error(err.into()),
            Some(Ok((chrom, size, mergingvalues))) => {
                let chrom_id = self.chrom_ids.as_mut().unwrap().get_id(&chrom);
                let group = BigWigWrite::begin_processing_chrom(
                    chrom,
                    chrom_id,
                    size,
                    mergingvalues,
                    self.pool.clone(),
                    self.options,
                );
                match group {
                    Ok(group) => ChromDataState::Read(group, self),
                    Err(err) => ChromDataState::Error(err.into()),
                }
            }
            None => {
                let chrom_ids = self.chrom_ids.take().unwrap();
                ChromDataState::Finished(chrom_ids)
            }
        }
    }
}

fn main() -> Result<(), WriteGroupsError> {
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
            let all_values = ChromGroupReadImpl {
                pool: futures::executor::ThreadPoolBuilder::new()
                    .pool_size(nthreads)
                    .create()
                    .expect("Unable to create thread pool."),
                options: outb.options,
                iter: Box::new(iter),
                chrom_ids: Some(IdMap::default()),
            };
            outb.write(chrom_map, all_values)?;
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
