#![feature(async_await, await_macro, futures_api, test)]

use std::io::{BufReader, BufRead};
use std::fs::File;

mod bigwig;
use bigwig::BigWigWrite;
use bigwig::ValueWithChrom;

mod idmap;
mod tell;

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
    println!("Path: {:?}", outb.path);

    let chrom_map = get_chrom_map(File::open(chroms)?);
    //println!("Chrom_map: {:?}", chrom_map);
    /*let chrom_map: std::collections::HashMap<&str, u32> = [
        ("chr2", 242193529),
        ("chr10", 133797422),
        ("chr17", 83257441),
        ("chr17_GL000205v2_random", 185591),
        ("chr17_KI270729v1_random", 280839),
        ("chr17_KI270730v1_random", 112551),
    ].iter().cloned().collect();*/

    println!("Reading file.");
    let infile = File::open(bg)?;
    //let infile = File::open("/home/hueyj/temp/test.bedGraph")?;
    let vals_iter = BufReader::new(infile)
        .lines()
        .map(|l| {
            let words = l.expect("Split error");
            let mut split = words.split_whitespace();
            ValueWithChrom {
                chrom: split.next().expect("Missing chrom").to_owned(),
                start: split.next().expect("Missing start").parse::<u32>().unwrap(),
                end: split.next().expect("Missing end").parse::<u32>().unwrap(),
                value: split.next().expect("Missing value").parse::<f32>().unwrap(),
            }
        });

    outb.write(chrom_map, vals_iter)?;
    Ok(())
}
