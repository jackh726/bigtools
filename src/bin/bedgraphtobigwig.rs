#![feature(async_await)]

use std::fs::File;
use std::io::{self, BufRead, BufReader};

use bigwig2::bigwig::BigWigWrite;
use bigwig2::bedgraphparser::{self, BedGraphParser};

fn main() -> io::Result<()> {
    let mut args = std::env::args();
    args.next();
    let bedgraphpath = args.next().expect("Must pass a bedgraph.");
    let chrom_map = args.next().expect("Must pass a chromosome sizes file.");
    let bigwigpath = args.next().expect("Must pass a bigwig.");
    println!("Args: {} {} {}", bedgraphpath, chrom_map, bigwigpath);


    let outb = BigWigWrite::create_file(bigwigpath)?;
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
