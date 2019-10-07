use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use clap::{App, Arg};

use bigwig2::bigwig::{BigWigWrite, WriteGroupsError};
use bigwig2::bedlikeparser;
use bigwig2::bedgraphparser::BedGraphParser;

fn main() -> Result<(), WriteGroupsError> {
    let matches = App::new("BedGraphToBigWig")
        .arg(Arg::with_name("bedgraph")
                .help("the bedgraph to convert to a bigwig")
                .index(1)
                .required(true)
            )
        .arg(Arg::with_name("chromsizes")
                .help("A chromosome sizes file. Each line should be have a chromosome and its size in bases, separated by whitespace.")
                .index(2)
                .required(true)
            )
        .arg(Arg::with_name("output")
                .help("The output bigwig path")
                .index(3)
                .required(true)
            )
        .get_matches();

    let bedgraphpath = matches.value_of("bedgraph").unwrap().to_owned();
    let chrom_map = matches.value_of("chromsizes").unwrap().to_owned();
    let bigwigpath = matches.value_of("output").unwrap().to_owned();

    let outb = BigWigWrite::create_file(bigwigpath);
    let chrom_map: HashMap<String, u32> = BufReader::new(File::open(chrom_map)?)
        .lines()
        .filter(|l| match l { Ok(s) => !s.is_empty(), _ => true })
        .map(|l| {
            let words = l.expect("Split error");
            let mut split = words.split_whitespace();
            (split.next().expect("Missing chrom").to_owned(), split.next().expect("Missing size").parse::<u32>().unwrap())
        })
        .collect();

    let pool = futures::executor::ThreadPoolBuilder::new().pool_size(6).create().expect("Unable to create thread pool.");

    let infile = File::open(bedgraphpath.clone())?;
    let vals_iter = BedGraphParser::from_file(infile);
    let options = outb.options.clone();

    let parse_fn = move |chrom, chrom_id, chrom_length, group| {
        BigWigWrite::begin_processing_chrom(chrom, chrom_id, chrom_length, group, pool.clone(), options.clone())
    };
    let chsi = bedlikeparser::BedGraphParserChromGroupStreamingIterator::new(vals_iter, chrom_map.clone(), Box::new(parse_fn));
    outb.write_groups(chrom_map, chsi)?;

    Ok(())
}
