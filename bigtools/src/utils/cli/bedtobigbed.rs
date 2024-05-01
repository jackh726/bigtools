use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use clap::Parser;
use tokio::runtime;

use crate::bed::bedparser::parse_bed;
use crate::bed::indexer::index_chroms;
use crate::bedchromdata::BedParserParallelStreamingIterator;
use crate::{
    bed::bedparser::BedParser, bedchromdata::BedParserStreamingIterator, BigBedWrite, InputSortType,
};

use super::BBIWriteArgs;

#[derive(Clone, Debug, PartialEq, Parser)]
#[command(
    name = "bedtobigbed",
    about = "Converts a bed to a bigBed.",
    long_about = None,
)]
pub struct BedToBigBedArgs {
    /// The bedGraph to convert to a bigbed.
    pub bed: String,

    /// A chromosome sizes file. Each line should be have a chromosome and its size in bases, separated by whitespace.
    pub chromsizes: String,

    /// The output bigwig path
    pub output: String,

    /// Path to a file containing the custom autosql to add to the bigBed file. If not specified, the standard BED
    /// autosql will be added (based on the number of columns).
    #[arg(short = 'a', long)]
    pub autosql: Option<String>,

    /// Set whether to read and convert the bedGraph in parallel. Requires that the bedGraph is sorted.
    /// Can take `auto` (default), `yes`, `no`. Ignored when input is stdin or when nthreads is `1`.
    #[arg(short = 'p', long)]
    #[arg(default_value = "auto")]
    pub parallel: String,

    /// If set, indicates that only a single pass should be done on the input file. This is most useful
    /// on large files in order to reduce total time. This automatically happens when the input is `stdin`.
    #[arg(long)]
    #[arg(default_value_t = false)]
    pub single_pass: bool,

    #[command(flatten)]
    pub write_args: BBIWriteArgs,
}

pub fn bedtobigbed(args: BedToBigBedArgs) -> Result<(), Box<dyn Error>> {
    let bedpath = args.bed;
    let chrom_map = args.chromsizes;
    let bigwigpath = args.output;
    let nthreads = args.write_args.nthreads;
    let input_sort_type = match args.write_args.sorted.as_ref() {
        "all" => InputSortType::ALL,
        "start" => InputSortType::START,
        "none" => {
            eprintln!("Using completely unsorted input is not implemented yet.");
            return Ok(());
        }
        sorted => {
            eprintln!(
                "Invalid option for `sorted`: `{}`. Options are `all`, `start`, or `none`.",
                sorted
            );
            return Ok(());
        }
    };

    let mut outb = BigBedWrite::create_file(bigwigpath);
    outb.options.max_zooms = args.write_args.nzooms;
    outb.options.compress = !args.write_args.uncompressed;
    outb.options.input_sort_type = input_sort_type;
    outb.options.inmemory = args.write_args.inmemory;
    let chrom_map: HashMap<String, u32> = BufReader::new(File::open(chrom_map)?)
        .lines()
        .filter(|l| match l {
            Ok(s) => !s.is_empty(),
            _ => true,
        })
        .map(|l| {
            let words = l.expect("Split error");
            let mut split = words.split_whitespace();
            (
                split.next().expect("Missing chrom").to_owned(),
                split.next().expect("Missing size").parse::<u32>().unwrap(),
            )
        })
        .collect();

    let runtime = if nthreads == 1 {
        outb.options.channel_size = 0;
        runtime::Builder::new_current_thread().build().unwrap()
    } else {
        runtime::Builder::new_multi_thread()
            .worker_threads(nthreads)
            .build()
            .unwrap()
    };

    let autosql = match args.autosql.as_ref() {
        None => {
            use crate::utils::chromvalues::ChromValues;
            let infile = File::open(&bedpath)?;
            let mut vals_iter = BedParser::from_bed_file(infile);
            let (_, mut group) = vals_iter.next_chrom().unwrap().unwrap();
            let first = group.peek().unwrap().unwrap();
            crate::bed::autosql::bed_autosql(&first.rest)
        }
        Some(file) => std::fs::read_to_string(file)?,
    };
    outb.autosql = Some(autosql);

    let allow_out_of_order_chroms = !matches!(outb.options.input_sort_type, InputSortType::ALL);
    if bedpath == "-" || bedpath == "stdin" {
        let stdin = std::io::stdin().lock();
        let vals_iter = BedParser::from_bed_file(stdin);

        let chsi = BedParserStreamingIterator::new(vals_iter, allow_out_of_order_chroms);
        outb.write_singlethreaded(chrom_map, chsi, runtime)?;
    } else {
        let infile = File::open(&bedpath)?;
        let (parallel, parallel_required) = match (nthreads, args.parallel.as_ref()) {
            (1, _) | (_, "no") => (false, false),
            (_, "auto") => (infile.metadata()?.len() >= 200_000_000, false),
            (_, "yes") => (true, true),
            (_, v) => {
                eprintln!(
                    "Unexpected value for `parallel`: \"{}\". Defaulting to `auto`.",
                    v
                );
                (infile.metadata()?.len() >= 200_000_000, false)
            }
        };
        let chrom_indices = match parallel {
            false => None,
            true => {
                let index = index_chroms(infile)?;
                match (index, parallel_required) {
                    (Some(index), _) => Some(index),
                    (None, true) => {
                        eprintln!(
                            "Parallel conversion requires a sorted bedGraph file. Cancelling.",
                        );
                        return Ok(());
                    }
                    (None, false) => None,
                }
            }
        };
        if let Some(chrom_indices) = chrom_indices {
            if args.single_pass {
                let chsi = BedParserParallelStreamingIterator::new(
                    chrom_indices,
                    allow_out_of_order_chroms,
                    PathBuf::from(bedpath),
                    parse_bed,
                );
                outb.write(chrom_map, chsi, runtime)?;
            } else {
                outb.write_multipass(
                    || {
                        let chsi = BedParserParallelStreamingIterator::new(
                            chrom_indices.clone(),
                            allow_out_of_order_chroms,
                            PathBuf::from(bedpath.clone()),
                            parse_bed,
                        );

                        Ok(chsi)
                    },
                    chrom_map,
                    runtime,
                )?;
            }
        } else {
            let infile = File::open(&bedpath)?;
            if args.single_pass {
                let vals_iter = BedParser::from_bed_file(infile);

                let chsi = BedParserStreamingIterator::new(vals_iter, allow_out_of_order_chroms);
                outb.write(chrom_map, chsi, runtime)?;
            } else {
                outb.write_multipass(
                    || {
                        let infile = File::open(&bedpath)?;
                        let vals_iter = BedParser::from_bed_file(infile);
                        let chsi =
                            BedParserStreamingIterator::new(vals_iter, allow_out_of_order_chroms);

                        Ok(chsi)
                    },
                    chrom_map,
                    runtime,
                )?;
            }
        }
    };

    Ok(())
}
