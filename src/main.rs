extern crate std;

mod bigwig;
use bigwig::BigWig;
use bigwig::ValueWithChrom;

fn main() -> Result<(), std::io::Error> {
    //write_test()?;

    //let mut b = BigWig::from_file(String::from("/home/hueyj/temp/ENCFF609KNT.bigWig"))?;
    let mut b = BigWig::from_file(String::from("/home/hueyj/temp/final.min.chr17.bigWig"))?;
    //let mut b = BigWig::from_file(String::from("/home/hueyj/temp/out.bigWig"))?;
    println!("Path: {:?}", b.path);

    b.read_info()?;

    //let interval = b.get_interval("chr1", 09000000u32, 10010000u32)?;
    //let interval = b.get_interval("chr17", 10000000u32, 10010000u32)?;
    //println!("Interval result: {:?} {:?}", interval.len(), &interval[0..10]);

    b.test_read_zoom("chr17", 10000000u32, 10010000u32)?;
    Ok(())
}

fn write_test() -> std::io::Result<()> {
    let mut outb = BigWig::create_file(String::from("/home/hueyj/temp/out.bigWig"))?;
    println!("Path: {:?}", outb.path);

    let chrom_map: std::collections::HashMap<&str, u32> = [
        ("chr2", 242193529),
        ("chr10", 133797422),
        ("chr17", 83257441)
    ].iter().cloned().collect();

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

    Ok(())
}
