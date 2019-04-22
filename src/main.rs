#![feature(async_await, await_macro, futures_api, test)]

use std::io::{BufReader, BufRead};
use std::fs::File;

mod bigwig;
use bigwig::BigWigRead;
use bigwig::BigWigWrite;

mod idmap;
mod tell;

mod bigwigmerge;
mod tempfilewrite;

fn main() -> Result<(), std::io::Error> {
    let mut args = std::env::args();
    args.next();
    let b1 = args.next().unwrap_or_else(|| "/home/hueyj/temp/final.min.chr17.bigWig".to_string());
    let b2 = args.next().unwrap_or_else(|| "/home/hueyj/temp/final.min.chr17.bigWig".to_string());
    println!("Args: {:} {:}", b1, b2);

    merge_test(b1, b2)?;

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

fn merge_test(b1path: String, b2path: String) -> std::io::Result<()> {
    //let b1 = BigWigRead::from_file_and_attach(String::from("/home/hueyj/temp/final.min.chr17.bigWig"))?;
    //let b2 = BigWigRead::from_file_and_attach(String::from("/home/hueyj/temp/final.min.chr17.bigWig"))?;
    //let mut b1 = BigWigRead::from_file_and_attach(String::from("/home/hueyj/temp/ENCFF207QBX.chr17.bigWig"))?;
    //let mut b2 = BigWigRead::from_file_and_attach(String::from("/home/hueyj/temp/ENCFF609KNT.chr17.bigWig"))?;

    let b1 = BigWigRead::from_file_and_attach(b1path)?;
    let b2 = BigWigRead::from_file_and_attach(b2path)?;

    let all_values = bigwigmerge::get_merged_values(vec![b1, b2])?;

    let out = String::from("/home/hueyj/temp/merge_test.bigWig");
    let outb = BigWigWrite::create_file(out)?;

    let chroms = String::from("/home/hueyj/temp/hg38.chrom.sizes");
    let chrom_map = get_chrom_map(File::open(chroms)?);

    outb.write(chrom_map, all_values)?;
    
    Ok(())
}
