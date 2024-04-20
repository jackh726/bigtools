use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

use clap::Parser;
use tokio::runtime;

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

    #[arg(short = 'a', long)]
    pub autosql: Option<String>,

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
        runtime::Builder::new_current_thread().build().unwrap()
    } else {
        runtime::Builder::new_multi_thread()
            .worker_threads(nthreads)
            .build()
            .unwrap()
    };

    let infile = File::open(bedpath)?;
    let mut vals_iter = BedParser::from_bed_file(infile);

    let autosql = match args.autosql.as_ref() {
        None => {
            let (_, mut group) = vals_iter.next_chrom().unwrap().unwrap();
            let first = group.peek_val().unwrap();
            crate::bed::autosql::bed_autosql(&first.rest)
        }
        Some(file) => std::fs::read_to_string(file)?,
    };
    outb.autosql = Some(autosql);

    let allow_out_of_order_chroms = !matches!(outb.options.input_sort_type, InputSortType::ALL);
    let chsi = BedParserStreamingIterator::new(vals_iter, allow_out_of_order_chroms);
    outb.write(chrom_map, chsi, runtime)?;

    Ok(())
}
