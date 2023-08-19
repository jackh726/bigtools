use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::FromStr;

use bigtools::bedchromdata::BedParserStreamingIterator;
use bigtools::utils::cli::BBIWriteArgs;
use clap::Parser;

use bigtools::bbi::BigBedWrite;
use bigtools::bbiwrite::InputSortType;
use bigtools::bed::bedparser::BedParser;

#[derive(Parser)]
#[command(about = "Converts a bed to a bigBed.", long_about = None)]
struct Cli {
    /// The n to convert to a bigbed.
    bed: String,

    /// A chromosome sizes file. Each line should be have a chromosome and its size in bases, separated by whitespace.
    chromsizes: String,

    /// The output bigwig path
    output: String,

    #[arg(short = 'a', long)]
    autosql: Option<String>,

    #[command(flatten)]
    write_args: BBIWriteArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args_os().map(|a| {
        match a.to_str() {
            Some("-unc") => return OsString::from_str("--uncompressed").unwrap(),
            Some("-tab") => return OsString::from_str("").unwrap(),
            Some(b) if b.starts_with("-blockSize=") => {
                return OsString::from_str(&format!(
                    "--block-size={}",
                    b.replace("-blockSize=", "")
                ))
                .unwrap()
            }
            Some(b) if b.starts_with("-itemsPerSlot=") => {
                return OsString::from_str(&format!(
                    "--items-per-slot={}",
                    b.replace("-itemsPerSlot=", "")
                ))
                .unwrap()
            }
            Some(b) if b.starts_with("-as=") => {
                return OsString::from_str(&format!("--autosql={}", b.replace("-as=", ""))).unwrap()
            }
            Some(b) if b.starts_with("-type=") => return OsString::from_str("").unwrap(),
            Some("-extraIndex")
            | Some("-sizesIs2Bit")
            | Some("-sizesIsChromAliasBb")
            | Some("-sizesIsBb")
            | Some("-allow1bpOverlap") => {
                panic!(
                    "Unimplemented compatibility option {}.",
                    a.to_string_lossy()
                );
            }
            Some(b) if b.starts_with("-extraIndex=") || b.starts_with("-udcDir") => {
                panic!(
                    "Unimplemented compatibility option {}.",
                    a.to_string_lossy()
                );
            }
            _ => {}
        }
        a
    });
    let matches = Cli::parse_from(args);

    let bedpath = matches.bed;
    let chrom_map = matches.chromsizes;
    let bigwigpath = matches.output;
    let nthreads = matches.write_args.nthreads;
    let input_sort_type = match matches.write_args.sorted.as_ref() {
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
    outb.options.max_zooms = matches.write_args.nzooms;
    outb.options.compress = !matches.write_args.uncompressed;
    outb.options.input_sort_type = input_sort_type;
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

    let pool = futures::executor::ThreadPoolBuilder::new()
        .pool_size(nthreads)
        .create()
        .expect("Unable to create thread pool.");

    let infile = File::open(bedpath)?;
    let mut vals_iter = BedParser::from_bed_file(infile);

    let autosql = match matches.autosql.as_ref() {
        None => {
            use bigtools::utils::chromvalues::ChromValues;
            let (_, mut group) = vals_iter.next_chrom().unwrap().unwrap();
            let first = group.peek().unwrap().unwrap();
            bigtools::bed::autosql::bed_autosql(&first.rest)
        }
        Some(file) => std::fs::read_to_string(file)?,
    };
    outb.autosql = Some(autosql);

    let allow_out_of_order_chroms = !matches!(outb.options.input_sort_type, InputSortType::ALL);
    let chsi = BedParserStreamingIterator::new(vals_iter, allow_out_of_order_chroms);
    outb.write(chrom_map, chsi, pool)?;

    Ok(())
}

#[test]
fn verify_cli_bedtobigbed() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
