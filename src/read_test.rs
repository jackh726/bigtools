#![feature(async_await, await_macro, test)]

mod bigwig;
use bigwig::BigWigRead;

mod idmap;
mod tell;
mod tempfilewrite;
mod bedgraphreader;
mod tempfilebuffer;

fn main() -> Result<(), std::io::Error> {
    let mut args = std::env::args();
    args.next();
    let out_file = args.next().unwrap_or_else(|| "/home/hueyj/temp/out.bigWig".to_string());
    println!("Args: {:}", out_file);

    read_test(out_file)?;

    Ok(())
}

fn read_test(bw: String) -> std::io::Result<()> {
    //let mut b = BigWig::from_file(String::from("/home/hueyj/temp/ENCFF609KNT.bigWig"))?;
    //let mut b = BigWig::from_file(String::from("/home/hueyj/temp/final.min.chr17.bigWig"))?;
    let b = BigWigRead::from_file_and_attach(bw)?;
    println!("Read path: {:?}", b.path);

    println!("BigWigInfo: {:?}", b.info);

    //let interval = b.get_interval("chr1", 09000000u32, 10010000u32)?;
    let interval = b.get_interval("chr17", 10_000_000u32, 10_010_000u32)?;
    //let interval = b.get_interval("chr17", 60000u32, 62000u32)?;
    println!("Interval result: {:?}", interval.collect::<Vec<_>>().len());

    b.test_read_zoom("chr17", 0, 83257441)?;
    //b.test_read_zoom("chr17", 10_000_000u32, 10_010_000u32)?;
    Ok(())
}
