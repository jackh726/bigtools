use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::str::FromStr;

use bigtools::bed::indexer::index_chroms;
use bigtools::bedchromdata::{BedParserParallelStreamingIterator, BedParserStreamingIterator};
use clap::Parser;

use bigtools::bbi::BigWigWrite;
use bigtools::bbiwrite::InputSortType;
use bigtools::bbiwrite::{DEFAULT_BLOCK_SIZE, DEFAULT_ITEMS_PER_SLOT};
use bigtools::bed::bedparser::{parse_bedgraph, BedParser};

#[derive(Debug, Parser)]
#[command(about = "Converts an input bedGraph to a bigWig. Can be multi-threaded for substantial speedups. Note that ~11 temporary files are created/maintained.", long_about = None)]
struct Cli {
    /// The bedgraph to convert to a bigwig. Can use `-` or `stdin` to read from stdin.
    bedgraph: String,

    /// A chromosome sizes file. Each line should be have a chromosome and its size in bases, separated by whitespace.
    chromsizes: String,

    /// The output bigwig path
    output: String,

    /// Set the number of threads to use. This tool will typically use ~225% CPU on a HDD. SDDs may be higher. (IO bound)
    #[arg(short = 't', long)]
    #[arg(default_value_t = 6)]
    nthreads: usize,

    /// Set the maximum of zooms to create.
    #[arg(short = 'z', long)]
    #[arg(default_value_t = 10)]
    nzooms: u32,

    /// Don't use compression.
    #[arg(short = 'u', long)]
    #[arg(default_value_t = false)]
    uncompressed: bool,

    /// Sets whether the input is sorted. Can take `all`, `start`, or `none`.
    /// `all` means that the input bedGraph is sorted by chroms and start (`sort -k1,1 -k2,2n`).
    /// `start` means that the the chroms are out of order but the starts within a chrom is sorted.
    /// `none` means that the file is not sorted at all.
    /// `all` is default. `none` currently errors but may be supported in the future.
    /// Note that using a value other than `all` will not guarantee (though likely) support for third-party tools.
    #[arg(short = 's', long)]
    #[arg(default_value = "all")]
    sorted: String,

    /// Set whether to read and convert the bedGraph in parallel.
    /// Can take `auto` (default), `yes`, `no`. Ignored when input is stdin or when nthreads is `1`.
    #[arg(short = 'p', long)]
    #[arg(default_value = "auto")]
    parallel: String,

    /// Number of items to bundle in r-tree.
    #[arg(long)]
    #[arg(default_value_t = DEFAULT_BLOCK_SIZE)]
    block_size: u32,

    /// Number of data points bundled at lowest level.
    #[arg(long)]
    #[arg(default_value_t = DEFAULT_ITEMS_PER_SLOT)]
    items_per_slot: u32,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args_os().map(|a| {
        match a.to_str() {
            Some("-unc") => return OsString::from_str("--uncompressed").unwrap(),
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
            _ => {}
        }
        a
    });
    let matches = Cli::parse_from(args);

    let bedgraphpath = matches.bedgraph;
    let chrom_map = matches.chromsizes;
    let bigwigpath = matches.output;
    let nthreads = matches.nthreads;
    let input_sort_type = match matches.sorted.as_ref() {
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

    let mut outb = BigWigWrite::create_file(bigwigpath);
    outb.options.max_zooms = matches.nzooms;
    outb.options.compress = !matches.uncompressed;
    outb.options.input_sort_type = input_sort_type;
    outb.options.block_size = matches.block_size;
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

    let allow_out_of_order_chroms = !matches!(outb.options.input_sort_type, InputSortType::ALL);
    if bedgraphpath == "-" || bedgraphpath == "stdin" {
        let stdin = std::io::stdin();
        // FIXME: This will lock on every line read, when we should be able to lock once
        let vals_iter = BedParser::from_bedgraph_file(stdin);

        let chsi = BedParserStreamingIterator::new(vals_iter, allow_out_of_order_chroms);
        outb.write(chrom_map, chsi, pool)?;
    } else {
        let infile = File::open(&bedgraphpath)?;
        let parallel = match (nthreads, matches.parallel.as_ref()) {
            (1, _) | (_, "auto") => infile.metadata()?.len() >= 200_000_000,
            (_, "yes") => true,
            (_, "no") => false,
            (_, v) => {
                eprintln!(
                    "Unexpected value for `parallel`: \"{}\". Defaulting to `auto`.",
                    v
                );
                true
            }
        };
        if parallel {
            let chrom_indices: Vec<(u64, String)> = index_chroms(infile)?;

            let chsi = BedParserParallelStreamingIterator::new(
                chrom_indices,
                allow_out_of_order_chroms,
                PathBuf::from(bedgraphpath),
                parse_bedgraph,
            );
            outb.write(chrom_map, chsi, pool)?;
        } else {
            let vals_iter = BedParser::from_bedgraph_file(infile);

            let chsi = BedParserStreamingIterator::new(vals_iter, allow_out_of_order_chroms);
            outb.write(chrom_map, chsi, pool)?;
        }
    };

    Ok(())
}

#[test]
fn verify_cli_bedgraphtobigwig() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
