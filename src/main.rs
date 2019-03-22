extern crate std;

use std::io::{BufReader, BufRead};
use std::fs::File;

mod bigwig;
use bigwig::BigWig;
use bigwig::ValueWithChrom;

mod idmap;

fn main() -> Result<(), std::io::Error> {
    write_test()?;

    //read_test()?;
    Ok(())
}

fn write_test() -> std::io::Result<()> {
    let mut outb = BigWig::create_file(String::from("/home/hueyj/temp/out.bigWig"))?;
    println!("Path: {:?}", outb.path);

    let chrom_map: std::collections::HashMap<&str, u32> = [
        ("chr2", 242193529),
        ("chr10", 133797422),
        ("chr17", 83257441),
        ("chr17_GL000205v2_random", 185591),
        ("chr17_KI270729v1_random", 280839),
        ("chr17_KI270730v1_random", 112551),
    ].iter().cloned().collect();

    println!("Reading file.");
    let infile = File::open("/home/hueyj/temp/final.min.chr17.full.bedGraph")?;
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

fn read_test() -> std::io::Result<()> {
    //let mut b = BigWig::from_file(String::from("/home/hueyj/temp/ENCFF609KNT.bigWig"))?;
    //let mut b = BigWig::from_file(String::from("/home/hueyj/temp/final.min.chr17.bigWig"))?;
    //let mut b = BigWig::from_file(String::from("/home/hueyj/temp/out.bigWig"))?;
    //println!("Path: {:?}", b.path);

    //b.read_info()?;

    //let interval = b.get_interval("chr1", 09000000u32, 10010000u32)?;
    //let interval = b.get_interval("chr17", 10000000u32, 10010000u32)?;
    //println!("Interval result: {:?} {:?}", interval.len(), &interval[0..10]);

    //b.test_read_zoom("chr17", 10000000u32, 10010000u32)?;
    Ok(())
}