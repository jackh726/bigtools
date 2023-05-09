use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use bigtools::bed::indexer::index_chroms;
use clap::{App, Arg};

use bigtools::bbiwrite::InputSortType;
use bigtools::bed::bedparser::{self, parse_bedgraph, BedParser};
use bigtools::bigwig::BigWigWrite;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("BedGraphToBigWig")
        .about("Converts an input bedGraph to a bigWig. Can be multi-threaded for substantial speedups. Note that ~11 temporary files are created/maintained.")
        .arg(Arg::new("bedgraph")
                .help("The bedgraph to convert to a bigwig. Can use `-` or `stdin` to read from stdin.")
                .index(1)
                .required(true)
            )
        .arg(Arg::new("chromsizes")
                .help("A chromosome sizes file. Each line should be have a chromosome and its size in bases, separated by whitespace.")
                .index(2)
                .required(true)
            )
        .arg(Arg::new("output")
                .help("The output bigwig path")
                .index(3)
                .required(true)
            )
        .arg(Arg::new("nthreads")
                .short('t')
                .help("Set the number of threads to use. This tool will typically use ~225% CPU on a HDD. SDDs may be higher. (IO bound)")
                .takes_value(true)
                .default_value("6"))
        .arg(Arg::new("nzooms")
                .short('z')
                .help("Set the maximum of zooms to create.")
                .takes_value(true)
                .default_value("10"))
        .arg(Arg::new("uncompressed")
                .short('u')
                .help("Don't use compression."))
        .arg(Arg::new("sorted")
                .short('s')
                .help("Sets whether the input is sorted. Can take `all`, `start`, or `none`. `all` means that the input bedGraph is sorted by chroms and start (`sort -k1,1 -k2,2n`). `start` means that the the chroms are out of order but the starts within a chrom is sorted. `none` means that the file is not sorted at all. `all` is default. `none` currently errors but may be supported in the future. Note that using a value other than `all` will not guarantee (though likely) support for third-party tools.")
                .takes_value(true)
                .default_value("all"))
        .get_matches();

    let bedgraphpath = matches.value_of("bedgraph").unwrap().to_owned();
    let chrom_map = matches.value_of("chromsizes").unwrap().to_owned();
    let bigwigpath = matches.value_of("output").unwrap().to_owned();
    let nthreads = {
        let nthreads = matches.value_of("nthreads").unwrap();
        match nthreads.parse() {
            Ok(parsed) => parsed,
            Err(_) => {
                eprintln!("Invalid argument for `nthreads`: must be a positive number");
                return Ok(());
            }
        }
    };
    let nzooms = {
        let nzooms = matches.value_of("nzooms").unwrap();
        match nzooms.parse() {
            Ok(parsed) => parsed,
            Err(_) => {
                eprintln!("Invalid argument for `nzooms`: must be a positive number");
                return Ok(());
            }
        }
    };
    let uncompressed = matches.is_present("uncompressed");
    let input_sort_type = match matches.value_of("sorted") {
        None => InputSortType::ALL,
        Some("all") => InputSortType::ALL,
        Some("start") => InputSortType::START,
        Some("none") => {
            eprintln!("Using completely unsorted input is not implemented yet.");
            return Ok(());
        }
        Some(sorted) => {
            eprintln!(
                "Invalid option for `sorted`: `{}`. Options are `all`, `start`, or `none`.",
                sorted
            );
            return Ok(());
        }
    };

    let mut outb = BigWigWrite::create_file(bigwigpath);
    outb.options.max_zooms = nzooms;
    outb.options.compress = !uncompressed;
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

    let allow_out_of_order_chroms = !matches!(outb.options.input_sort_type, InputSortType::ALL);
    if bedgraphpath == "-" || bedgraphpath == "stdin" {
        let stdin = std::io::stdin();
        // FIXME: This will lock on every line read, when we should be able to lock once
        let vals_iter = BedParser::from_bedgraph_file(stdin);

        let chsi = bedparser::BedParserStreamingIterator::new(
            vals_iter,
            chrom_map.clone(),
            allow_out_of_order_chroms,
        );
        outb.write(chrom_map, chsi, pool)?;
    } else {
        let infile = File::open(&bedgraphpath)?;
        if infile.metadata()?.len() >= 200_000_000 {
            let chrom_indices: Vec<(u64, String)> = index_chroms(infile)?;

            let chsi = bedparser::BedParserParallelStreamingIterator::new(
                chrom_map.clone(),
                chrom_indices,
                allow_out_of_order_chroms,
                PathBuf::from(bedgraphpath),
                parse_bedgraph,
            );
            outb.write(chrom_map, chsi, pool)?;
        } else {
            let vals_iter = BedParser::from_bedgraph_file(infile);

            let chsi = bedparser::BedParserStreamingIterator::new(
                vals_iter,
                chrom_map.clone(),
                allow_out_of_order_chroms,
            );
            outb.write(chrom_map, chsi, pool)?;
        }
    };

    Ok(())
}
