use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use clap::Parser;
use crossbeam_channel::unbounded;
use thiserror::Error;

use crate::utils::chromvalues::ChromValues;
use crate::utils::merge::merge_sections_many;
use crate::utils::reopen::ReopenableFile;
use crate::Value;
use crate::{BBIReadError, BigWigRead, BigWigWrite};
use crate::{ChromData, ChromDataState, ChromProcessingKey, ProcessChromError};
use tokio::runtime;

use super::BBIWriteArgs;

#[derive(Clone, Debug, PartialEq, Parser)]
#[command(
    name = "bigwigmerge",
    about = "Merges multiple bigwigs.",
    long_about = None,
)]
pub struct BigWigMergeArgs {
    /// the path of the merged output bigwig (if .bw or .bigWig) or bedGraph (if .bedGraph)
    pub output: String,

    /// the path of an input bigwig to merge
    #[arg(short = 'b')]
    pub bigwig: Vec<String>,

    /// a line-delimited list of bigwigs
    #[arg(short = 'l')]
    pub list: Vec<String>,

    /// Don't output values at or below this threshold. Default is 0.0
    #[arg(long)]
    #[arg(default_value_t = 0.0)]
    pub threshold: f32,

    /// Add adjustment to each value
    #[arg(long)]
    pub adjust: Option<f32>,

    /// Values higher than this are clipped to this value
    #[arg(long)]
    pub clip: Option<f32>,

    /// Merged value is maximum from input files rather than sum
    #[arg(long)]
    #[arg(default_value_t = false)]
    max: bool,

    /// Can be `bigwig` or `bedgraph` (case-insensitive). If not specified,
    /// will be inferred from the output file ending.
    #[arg(long)]
    output_type: Option<String>,

    #[command(flatten)]
    write_args: BBIWriteArgs,
}

pub fn bigwigmerge(args: BigWigMergeArgs) -> Result<(), Box<dyn Error>> {
    let output = args.output;
    let mut bigwigs: Vec<BigWigRead<ReopenableFile>> = vec![];

    for name in args.bigwig {
        match BigWigRead::open_file(&name) {
            Ok(bw) => bigwigs.push(bw),
            Err(e) => {
                eprintln!("Error when opening bigwig ({}): {:?}", name, e);
                return Ok(());
            }
        }
    }
    for list in args.list {
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

    let nthreads = args.write_args.nthreads;

    let (iter, chrom_map) = get_merged_vals(bigwigs, 10, args.threshold, args.adjust, args.clip)?;

    enum OutputType {
        BigWig,
        BedGraph,
    }

    let output_type = match (args.output_type, &output) {
        (None, output)
            if output.to_lowercase().ends_with(".bw")
                || output.to_lowercase().ends_with(".bigWig") =>
        {
            OutputType::BigWig
        }
        (None, output) if output.to_lowercase().ends_with(".bedGraph") => OutputType::BedGraph,
        (Some(output_type), _) if output_type.to_lowercase() == "bigwig" => OutputType::BigWig,
        (Some(output_type), _) if output_type.to_lowercase() == "bedgraph" => OutputType::BedGraph,
        _ => {
            eprintln!("Unable to determine output file format. \
                The output file must either in with `.bw` or `.bigWig` for bigwigs or `.bedGraph` for bedGraphs; or \
                `--output-type` must be set to either `bigwig` or `bedgraph`.");
            return Ok(());
        }
    };
    match output_type {
        OutputType::BigWig => {
            let outb = BigWigWrite::create_file(output);
            let runtime = if nthreads == 1 {
                runtime::Builder::new_current_thread().build().unwrap()
            } else {
                runtime::Builder::new_multi_thread()
                    .worker_threads(nthreads)
                    .build()
                    .unwrap()
            };
            let all_values = ChromGroupReadImpl {
                iter: Box::new(iter),
            };
            outb.write(chrom_map, all_values, runtime)?;
        }
        OutputType::BedGraph => {
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
    }

    //TODO: fails with too many open files
    Ok(())
}

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
        for chrom in bigwigs
            .iter()
            .flat_map(BigWigRead::chroms)
            .map(|c| c.name.clone())
        {
            if chrom_sizes.get(&chrom).is_some() {
                continue;
            }
            let mut size = None;
            let mut bws = Vec::with_capacity(bigwigs.len());
            for w in bigwigs.iter() {
                let chroms = w.chroms();
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
                bws.push((w.info().clone(), w.inner_read().path.to_string()));
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
