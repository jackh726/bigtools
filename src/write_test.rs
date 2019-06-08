#![feature(async_await, await_macro, test)]

use std::io::{BufReader, BufRead};
use std::fs::File;

mod bigwig;
use bigwig::BigWigWrite;
use bigwig::ValueWithChrom;

mod idmap;
mod tell;
mod tempfilewrite;
mod bedgraphreader;
mod tempfilebuffer;
mod streaming_linereader;
mod bedgraphparser;
mod chromvalues;

fn main() -> Result<(), std::io::Error> {
    let mut args = std::env::args();
    args.next();
    let bedgraph = args.next().unwrap_or_else(|| "/home/hueyj/temp/final.min.chr17.full.bedGraph".to_string());
    let chrom_map = args.next().unwrap_or_else(|| "/home/hueyj/temp/hg38.chrom.sizes".to_string());
    let out_file = args.next().unwrap_or_else(|| "/home/hueyj/temp/out.bigWig".to_string());
    println!("Args: {:} {:} {:}", bedgraph, chrom_map, out_file);

    write_test(bedgraph, chrom_map, out_file.clone())?;

    Ok(())
}

fn get_chrom_map(file: File) -> std::collections::HashMap<String, u32> {
    BufReader::new(file)
        .lines()
        .filter(|l| match l { Ok(s) => !s.is_empty(), _ => true })
        .map(|l| {
            let words = l.expect("Split error");
            let mut split = words.split_whitespace();
            (split.next().expect("Missing chrom").to_owned(), split.next().expect("Missing size").parse::<u32>().unwrap())
        })
        .collect()
}

fn write_test(bg: String, chroms: String, out: String) -> std::io::Result<()> {
    let outb = BigWigWrite::create_file(out)?;

    let chrom_map = get_chrom_map(File::open(chroms)?);

    println!("Reading file.");
    let infile = File::open(bg)?;

    let vals_iter = bedgraphparser::BedGraphParser::new(infile);
    outb.write(chrom_map, vals_iter)?;
    Ok(())
}
