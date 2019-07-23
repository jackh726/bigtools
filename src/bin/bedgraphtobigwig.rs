#![feature(async_await)]

use std::fs::File;
use std::io::{BufRead, BufReader};

use clap::{App, Arg};

use bigwig2::bigwig::{BigWigWrite, WriteGroupsError};
use bigwig2::bedgraphparser::{self, BedGraphParser};

fn main() -> Result<(), WriteGroupsError> {
    let matches = App::new("BigWigInfo")
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
    let chrom_map = BufReader::new(File::open(chrom_map)?)
        .lines()
        .filter(|l| match l { Ok(s) => !s.is_empty(), _ => true })
        .map(|l| {
            let words = l.expect("Split error");
            let mut split = words.split_whitespace();
            (split.next().expect("Missing chrom").to_owned(), split.next().expect("Missing size").parse::<u32>().unwrap())
        })
        .collect();

    let infile = File::open(bedgraphpath)?;
    let vals_iter = BedGraphParser::<BufReader<File>>::new(infile);
    let chsi = bedgraphparser::get_chromgroupstreamingiterator(vals_iter, outb.options.clone());
    outb.write_groups(chrom_map, chsi)?;

    Ok(())
}
