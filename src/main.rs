extern crate std;

use std::io::{BufReader, BufRead};
use std::fs::File;

mod bigwig;
use bigwig::BigWig;
use bigwig::ValueWithChrom;

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
/*
    let mut vals: Vec<ValueWithChrom> = vec![];
    for i in 0..500000 {
        vals.push(ValueWithChrom {
            chrom: "chr17",
            start: i + 10000000u32,
            end: i + 10000000u32 + 1,
            value: (i as f32 / 500000 as f32),
        });
    }
    vals.push(ValueWithChrom {
        chrom: "chr10",
        start: 1,
        end: 50,
        value: 0.3,
    });
    vals.push(ValueWithChrom {
        chrom: "chr10",
        start: 51,
        end: 100,
        value: 0.4,
    });
    // for i in 1..=10000 {
    //     vals.push(ValueWithChrom {
    //         chrom: "chr2",
    //         start: i,
    //         end: i+1,
    //         value: (i as f32 / 10000 as f32),
    //     });
    // }
    vals.sort_by_key(|v| (v.chrom, v.start));
    outb.write(chrom_map, vals.into_iter())?;
*/

    //println!("Reading line count");
    //let mut infile = File::open("/home/hueyj/temp/final.min.chr17.bedGraph")?;
    //let numlines = {
    //    let f = BufReader::new(infile);
    //    let numlines = f.lines().count();
    //    numlines
    //};
    let numlines = 78722367;
    //let numlines = 1000000;

    println!("Reading file.");
    let mut infile = File::open("/home/hueyj/temp/final.min.chr17.bedGraph")?;
    //let infile = File::open("/home/hueyj/temp/test.bedGraph")?;
    let mut vals_iter = BufReader::new(infile)
        .lines()
        .map(|l| {
            let words = l.expect("Split error");
            let split: Vec<String> = words.split_whitespace().map(|w| w.to_owned()).collect();
            ValueWithChrom {
                chrom: split[0].clone(),
                start: split[1].clone().parse::<u32>().unwrap(),
                end: split[2].clone().parse::<u32>().unwrap(),
                value: split[3].clone().parse::<f32>().unwrap(),
            }
        });

    use indicatif::ProgressBar;
    let pb = ProgressBar::new(numlines as u64);
    pb.set_draw_delta(50000);
    let progress_iter = pb.wrap_iter(vals_iter);

    outb.write(chrom_map, progress_iter)?;
    Ok(())
}

fn read_test() -> std::io::Result<()> {
    //let mut b = BigWig::from_file(String::from("/home/hueyj/temp/ENCFF609KNT.bigWig"))?;
    //let mut b = BigWig::from_file(String::from("/home/hueyj/temp/final.min.chr17.bigWig"))?;
    let mut b = BigWig::from_file(String::from("/home/hueyj/temp/out.bigWig"))?;
    //println!("Path: {:?}", b.path);

    b.read_info()?;

    //let interval = b.get_interval("chr1", 09000000u32, 10010000u32)?;
    //let interval = b.get_interval("chr17", 10000000u32, 10010000u32)?;
    //println!("Interval result: {:?} {:?}", interval.len(), &interval[0..10]);

    //b.test_read_zoom("chr17", 10000000u32, 10010000u32)?;
    Ok(())
}