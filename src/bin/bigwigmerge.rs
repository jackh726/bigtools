use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use clap::{App, Arg};

use futures::future::Either;

use bigtools::bigwig::BBIWriteOptions;
use bigtools::bigwig::ChromGroupRead;
use bigtools::bigwig::ChromGroupReadStreamingIterator;
use bigtools::bigwig::Value;
use bigtools::bigwig::{BBIRead, BigWigRead, BigWigWrite, WriteGroupsError};
use bigtools::chromvalues::ChromValues;

use bigtools::idmap::IdMap;
use bigtools::seekableread::ReopenableFile;
use bigtools::utils::merge_sections_many;

pub struct MergingValues {
    iter: std::iter::Peekable<Box<dyn Iterator<Item = Value> + Send>>,
}

impl ChromValues<Value> for MergingValues {
    fn next(&mut self) -> io::Result<Option<Value>> {
        Ok(self.iter.next())
    }

    fn peek(&mut self) -> Option<&Value> {
        self.iter.peek()
    }
}

pub fn get_merged_vals(
    bigwigs: Vec<BigWigRead<ReopenableFile, File>>,
) -> io::Result<(
    impl Iterator<Item = io::Result<(String, u32, MergingValues)>>,
    HashMap<String, u32>,
)> {
    // Get sizes for each and check that all files (that have the chrom) agree
    // Check that all chrom sizes match for all files
    let mut chrom_sizes = BTreeMap::new();
    let mut chrom_map = HashMap::new();
    for chrom in bigwigs.iter().flat_map(BBIRead::get_chroms).map(|c| c.name) {
        if chrom_sizes.get(&chrom).is_some() {
            continue;
        }
        let (sizes, bws): (Vec<_>, Vec<_>) = bigwigs
            .iter()
            .map(|w| {
                let chroms = w.get_chroms();
                let res = chroms.iter().find(|v| v.name == chrom);
                match res {
                    Some(s) => Some((s.length, w.clone())),
                    None => None,
                }
            })
            .filter_map(|x| x)
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

    let iter = chrom_sizes.clone().into_iter().map(|(chrom, (size, bws))| {
        let iters: Vec<_> = bws
            .into_iter()
            .map(|b| b.get_interval_move(&chrom, 1, size))
            .collect::<io::Result<Vec<_>>>()?;
        let iter: Box<dyn Iterator<Item = Value> + Send> =
            Box::new(merge_sections_many(iters).filter(|x| x.value != 0.0));
        let mergingvalues = MergingValues {
            iter: iter.peekable(),
        };

        Ok((chrom, size, mergingvalues))
    });

    Ok((iter, chrom_map))
}

pub fn get_merged_values(
    iter: impl Iterator<Item = io::Result<(String, u32, MergingValues)>> + Send + 'static,
    options: BBIWriteOptions,
    nthreads: usize,
) -> io::Result<impl ChromGroupReadStreamingIterator + std::marker::Send> {
    struct ChromGroupReadStreamingIteratorImpl {
        pool: futures::executor::ThreadPool,
        options: BBIWriteOptions,
        iter: Box<dyn Iterator<Item = io::Result<(String, u32, MergingValues)>> + Send>,
        chrom_ids: Option<IdMap>,
    }

    impl ChromGroupReadStreamingIterator for ChromGroupReadStreamingIteratorImpl {
        fn next(&mut self) -> Result<Option<Either<ChromGroupRead, (IdMap)>>, WriteGroupsError> {
            let next = self.iter.next();
            match next {
                Some(next) => {
                    let (chrom, size, mergingvalues) = next?;
                    let chrom_id = self.chrom_ids.as_mut().unwrap().get_id(&chrom);
                    let group = BigWigWrite::begin_processing_chrom(
                        chrom,
                        chrom_id,
                        size,
                        mergingvalues,
                        self.pool.clone(),
                        self.options.clone(),
                    )?;
                    Ok(Some(Either::Left(group)))
                }
                None => match self.chrom_ids.take() {
                    Some(chrom_ids) => Ok(Some(Either::Right(chrom_ids))),
                    None => Ok(None),
                },
            }
        }
    }

    let group_iter = ChromGroupReadStreamingIteratorImpl {
        pool: futures::executor::ThreadPoolBuilder::new()
            .pool_size(nthreads)
            .create()
            .expect("Unable to create thread pool."),
        options,
        iter: Box::new(iter),
        chrom_ids: Some(IdMap::default()),
    };

    Ok(group_iter)
}

fn main() -> Result<(), WriteGroupsError> {
    let matches = App::new("BigWigMerge")
        .arg(Arg::with_name("output")
                .help("the path of the merged output bigwig (if .bw or .bigWig) or bedGraph (if .bedGraph)")
                .index(1)
                .required(true)
            )
        .arg(Arg::with_name("bigwig")
                .short("b")
                .help("the path of an input bigwig to merge")
                .multiple(true)
                .takes_value(true)
            )
        .arg(Arg::with_name("list")
                .short("l")
                .help("a line-delimited list of bigwigs")
                .multiple(true)
                .takes_value(true)
            )
        .arg(Arg::with_name("nthreads")
                .short("t")
                .help("Set the number of threads to use")
                .takes_value(true)
                .default_value("6"))
        .get_matches();

    let output = matches.value_of("output").unwrap().to_owned();
    let mut bigwigs: Vec<BigWigRead<ReopenableFile, File>> = vec![];

    if let Some(bws) = matches.values_of("bigwig") {
        let results = bws
            .map(|b| BigWigRead::from_file_and_attach(b.to_owned()))
            .collect::<Result<Vec<_>, _>>();
        match results {
            Ok(bws) => bigwigs.extend(bws),
            Err(e) => {
                eprintln!("Error: {:?}", e);
                return Ok(());
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
                match BigWigRead::from_file_and_attach(name.clone()) {
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

    let (iter, chrom_map) = get_merged_vals(bigwigs)?;

    match output {
        output if output.ends_with(".bw") || output.ends_with(".bigWig") => {
            let outb = BigWigWrite::create_file(output);
            let all_values = get_merged_values(iter, outb.options.clone(), nthreads)?;
            outb.write_groups(chrom_map, all_values)?;
        }
        output if output.ends_with(".bedGraph") => {
            // TODO: convert to multi-threaded
            use std::io::Write;

            let mut chroms: Vec<String> = chrom_map.keys().map(|c| c.to_string()).collect();
            chroms.sort();

            let bedgraph = File::create(output)?;
            let mut writer = io::BufWriter::new(bedgraph);

            for v in iter {
                let (chrom, _, mut values) = v?;
                while let Some(val) = values.next()? {
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
