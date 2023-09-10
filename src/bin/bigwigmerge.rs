use std::collections::{BTreeMap, HashMap};
use std::env;
use std::error::Error;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::str::FromStr;

use clap::Parser;
use crossbeam_channel::unbounded;
use thiserror::Error;

use bigtools::utils::chromvalues::ChromValues;
use bigtools::utils::merge::merge_sections_many;
use bigtools::utils::reopen::ReopenableFile;
use bigtools::Value;
use bigtools::{BBIRead, BBIReadError, BigWigRead, BigWigWrite};
use bigtools::{ChromData, ChromDataState, ChromProcessingKey, ProcessChromError};

pub struct MergingValues {
    // We Box<dyn Iterator> because other this would be a mess to try to type
    iter: std::iter::Peekable<Box<dyn Iterator<Item = Result<Value, MergingValuesError>> + Send>>,
}

impl MergingValues {
    pub fn new<I: 'static>(
        iters: Vec<I>,
        threshold: f32,
        adjust: Option<f32>,
        clip: Option<f32>,
    ) -> Self
    where
        I: Iterator<Item = Result<Value, MergingValuesError>> + Send,
    {
        let adjust = adjust.unwrap_or(0.0);
        let iter: Box<dyn Iterator<Item = Result<Value, MergingValuesError>> + Send> = Box::new(
            merge_sections_many(iters)
                .map(move |x| {
                    x.map(|mut v| {
                        if let Some(clip) = clip {
                            v.value = clip.min(v.value);
                        }
                        v.value += adjust;
                        v
                    })
                })
                .filter(move |x| x.as_ref().map_or(true, |v| v.value > threshold)),
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
    bigwigs: Vec<BigWigRead<ReopenableFile>>,
    max_zooms: usize,
    threshold: f32,
    adjust: Option<f32>,
    clip: Option<f32>,
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
            let mut size = None;
            let mut bws = Vec::with_capacity(bigwigs.len());
            for w in bigwigs.iter() {
                let chroms = w.get_chroms();
                let res = chroms.iter().find(|v| v.name == chrom);
                let res = match res {
                    Some(res) => res,
                    None => continue,
                };
                match size {
                    Some(all_size) => {
                        if all_size != res.length {
                            eprintln!("Chrom '{:?}' had different sizes in the bigwig files. (Are you using the same assembly?)", chrom);
                            return Err(MergingValuesError::MismatchedChroms(
                                "Invalid input (nonmatching chroms)".to_owned(),
                            ));
                        }
                    }
                    None => {
                        size = Some(res.length);
                    }
                }
                // We don't want to a new file descriptor for every chrom
                bws.push((w.get_info().clone(), w.inner_read().path.to_string()));
            }
            let size = size.unwrap();

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
                    let f = ReopenableFile { file: File::open(&b.1)?, path: b.1 };
                    let b = BigWigRead::with_info(b.0, f);
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
                        let mut mergingvalues = MergingValues::new(chunk, threshold, adjust, clip);
                        let (sender, receiver) = unbounded::<Value>();
                        while let Some(val) = mergingvalues.next() {
                            let val = val?;
                            sender.send(val).unwrap();
                        }

                        merges.push(Box::new(receiver.into_iter().map(Result::Ok)));
                    }
                    merges
                };
            }

            let mergingvalues = MergingValues::new(merges, threshold, adjust, clip);
            Ok((chrom, size, mergingvalues))
        } else {
            let iters: Vec<_> = bws
                .into_iter()
                .map(|b| {
                    let f = ReopenableFile { file: File::open(&b.1)?, path: b.1 };
                    let b = BigWigRead::with_info(b.0, f);
                    b.get_interval_move(&chrom, 1, size).map(|i| i.map(|r| r.map_err(|e| MergingValuesError::BBIReadError(e))))
                })
                .collect::<Result<Vec<_>, _>>()?;
            let mergingvalues = MergingValues::new(iters, threshold, adjust, clip);

            Ok((chrom, size, mergingvalues))
        }
    });

    Ok((iter, chrom_map))
}

struct ChromGroupReadImpl {
    iter: Box<dyn Iterator<Item = Result<(String, u32, MergingValues), MergingValuesError>> + Send>,
}

impl ChromData for ChromGroupReadImpl {
    type Values = MergingValues;

    fn advance<
        State,
        F: FnMut(
            String,
            MergingValues,
            &mut State,
        ) -> Result<ChromProcessingKey, ProcessChromError<MergingValuesError>>,
    >(
        &mut self,
        do_read: &mut F,
        state: &mut State,
    ) -> Result<
        ChromDataState<ChromProcessingKey, MergingValuesError>,
        ProcessChromError<MergingValuesError>,
    > {
        let next: Option<Result<(String, u32, MergingValues), MergingValuesError>> =
            self.iter.next();
        Ok(match next {
            Some(Err(err)) => ChromDataState::Error(err.into()),
            Some(Ok((chrom, _, mergingvalues))) => {
                let read = do_read(chrom, mergingvalues, state)?;

                ChromDataState::NewChrom(read)
            }
            None => ChromDataState::Finished,
        })
    }
}

#[derive(Parser)]
#[command(about = "Merges multiple bigwigs.", long_about = None)]
struct Cli {
    /// the path of the merged output bigwig (if .bw or .bigWig) or bedGraph (if .bedGraph)
    output: String,

    /// the path of an input bigwig to merge
    #[arg(short = 'b')]
    bigwig: Vec<String>,

    /// a line-delimited list of bigwigs
    #[arg(short = 'l')]
    list: Vec<String>,

    /// Don't output values at or below this threshold. Default is 0.0
    #[arg(long)]
    #[arg(default_value_t = 0.0)]
    threshold: f32,

    /// Add adjustment to each value
    #[arg(long)]
    adjust: Option<f32>,

    /// Values higher than this are clipped to this value
    #[arg(long)]
    clip: Option<f32>,

    /// Merged value is maximum from input files rather than sum
    #[arg(long)]
    #[arg(default_value_t = false)]
    max: bool,

    /// Set the number of threads to use. This tool will nearly always benefit from more cores (<= # chroms).
    /// Note: for parts of the runtime, the actual usage may be nthreads+1
    #[arg(short = 't', long)]
    #[arg(default_value_t = 6)]
    nthreads: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = env::args_os().collect();
    let has_input = args.iter().any(|a| {
        a.to_str()
            .map_or(false, |a| a.starts_with("-b") || a.starts_with("-l"))
    });
    let args = if !has_input {
        // If there are no -l or -b, then let's see if it looks like a kent bigWigMerge call
        let in_list = args
            .iter()
            .any(|a| a.to_str().map_or(false, |a| a == "-inList"));
        let mut old_args = args;
        let mut args = vec![];
        old_args.reverse();
        while let Some(os_arg) = old_args.pop() {
            let arg = os_arg.to_string_lossy();
            if arg == "-inList" {
                continue;
            }
            if arg.starts_with("-") && !arg.contains("=") {
                args.push(os_arg);
                args.pop().map(|a| args.push(a));
                continue;
            }
            let more_args = 'more: {
                let mut args_iter = args.iter().rev().peekable();
                while let Some(arg) = args_iter.next() {
                    let arg = arg.to_string_lossy();
                    if arg.starts_with("-") && !arg.contains("=") {
                        args_iter.next();
                    } else {
                        break 'more true;
                    }
                }
                false
            };
            if !more_args {
                if in_list {
                    args.push(OsString::from_str("-l").unwrap());
                } else {
                    args.push(OsString::from_str("-b").unwrap());
                }
            }
            args.push(os_arg);
        }
        args
    } else {
        args
    };

    let args = args.into_iter().map(|a| {
        bigtools::compat_replace!(a;
            replace:
                "-threshold", "--threshold";
                "-adjust", "--adjust";
                "-clip", "--clip"
            ignore:
                "-inList"
            unimplemented:
                "-max";
                "-udcDir"
        )
    });
    let matches = Cli::parse_from(args);

    let output = matches.output;
    let mut bigwigs: Vec<BigWigRead<ReopenableFile>> = vec![];

    for name in matches.bigwig {
        match BigWigRead::open_file(&name) {
            Ok(bw) => bigwigs.push(bw),
            Err(e) => {
                eprintln!("Error when opening bigwig ({}): {:?}", name, e);
                return Ok(());
            }
        }
    }
    for list in matches.list {
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
            match BigWigRead::open_file(&name) {
                Ok(bw) => bigwigs.push(bw),
                Err(e) => {
                    eprintln!("Error when opening bigwig ({}): {:?}", name, e);
                    return Ok(());
                }
            }
        }
    }

    let nthreads = matches.nthreads;

    let (iter, chrom_map) =
        get_merged_vals(bigwigs, 10, matches.threshold, matches.adjust, matches.clip)?;

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

#[test]
fn verify_cli_bigwigmerge() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
