use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use bigtools::bed::indexer::index_chroms;
use bigtools::bedchromdata::{BedParserParallelStreamingIterator, BedParserStreamingIterator};
use bigtools::utils::cli::BBIWriteArgs;
use clap::Parser;

use bigtools::bbi::BigWigWrite;
use bigtools::bbiwrite::InputSortType;
use bigtools::bed::bedparser::{parse_bedgraph, BedParser};

#[derive(Parser)]
#[command(about = "Converts an input bedGraph to a bigWig. Can be multi-threaded for substantial speedups. Note that ~11 temporary files are created/maintained.", long_about = None)]
struct Cli {
    /// The bedgraph to convert to a bigwig. Can use `-` or `stdin` to read from stdin.
    bedgraph: String,

    /// A chromosome sizes file. Each line should be have a chromosome and its size in bases, separated by whitespace.
    chromsizes: String,

    /// The output bigwig path
    output: String,

    /// Set whether to read and convert the bedGraph in parallel.
    /// Can take `auto` (default), `yes`, `no`. Ignored when input is stdin or when nthreads is `1`.
    #[arg(short = 'p', long)]
    #[arg(default_value = "auto")]
    pub parallel: String,

    /// If set, indicates that only a single pass should be done on the input file. This is most useful
    /// on large files in order to reduce total time. This automatically happens when the input is `stdin`.
    #[arg(long)]
    #[arg(default_value_t = true)]
    single_pass: bool,

    #[command(flatten)]
    write_args: BBIWriteArgs,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args_os().map(|a| {
        bigtools::compat_replace!(a;
            replace:
                "-unc", "--uncompressed";
                "-blockSize", "--block-size";
                "-itemsPerSlot", "--items-per-slot"
            ignore:
            unimplemented:
        )
    });
    let matches = Cli::parse_from(args);

    let bedgraphpath = matches.bedgraph;
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

    let mut outb = BigWigWrite::create_file(bigwigpath);
    outb.options.max_zooms = matches.write_args.nzooms;
    outb.options.compress = !matches.write_args.uncompressed;
    outb.options.input_sort_type = input_sort_type;
    outb.options.block_size = matches.write_args.block_size;
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
        let stdin = std::io::stdin().lock();
        let vals_iter = BedParser::from_bedgraph_file(stdin);

        let chsi = BedParserStreamingIterator::new(vals_iter, allow_out_of_order_chroms);
        outb.write_singlethreaded(chrom_map, chsi, pool)?;
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

            if matches.single_pass {
                let chsi = BedParserParallelStreamingIterator::new(
                    chrom_indices,
                    allow_out_of_order_chroms,
                    PathBuf::from(bedgraphpath),
                    parse_bedgraph,
                );
                outb.write(chrom_map, chsi, pool)?;
            } else {
                outb.write_multipass(
                    || {
                        let chsi = BedParserParallelStreamingIterator::new(
                            chrom_indices.clone(),
                            allow_out_of_order_chroms,
                            PathBuf::from(bedgraphpath.clone()),
                            parse_bedgraph,
                        );

                        Ok(chsi)
                    },
                    chrom_map,
                    pool,
                )?;
            }
        } else {
            if matches.single_pass {
                let vals_iter = BedParser::from_bedgraph_file(infile);

                let chsi = BedParserStreamingIterator::new(vals_iter, allow_out_of_order_chroms);
                outb.write(chrom_map, chsi, pool)?;
            } else {
                outb.write_multipass(
                    || {
                        let infile = File::open(&bedgraphpath)?;
                        let vals_iter = BedParser::from_bedgraph_file(infile);
                        let chsi =
                            BedParserStreamingIterator::new(vals_iter, allow_out_of_order_chroms);

                        Ok(chsi)
                    },
                    chrom_map,
                    pool,
                )?;
            }
        }
    };

    Ok(())
}

#[test]
fn verify_cli_bedgraphtobigwig() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
