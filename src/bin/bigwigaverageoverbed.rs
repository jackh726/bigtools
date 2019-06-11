#![feature(async_await)]

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};

use bigwig2::bigwig::BigWigRead;
use bigwig2::streaming_linereader::StreamingLineReader;

fn write(bedin: File, bigwigin: BigWigRead, bedout: File) -> io::Result<()> {
    let mut bedstream = StreamingLineReader::new(BufReader::new(bedin));
    let mut bedoutwriter = BufWriter::new(bedout);

    while let Some(line) = bedstream.read()? {
        let mut split = line.splitn(4, "\t");
        let chrom = split.next().expect("Missing chrom");
        let start = split.next().expect("Missing start").parse::<u32>().unwrap();
        let end = split.next().expect("Missing end").parse::<u32>().unwrap();
        let rest = split.next().unwrap_or("").trim_end();
        let interval = bigwigin.get_interval(chrom, start, end)?.collect::<Result<Vec<_>, _>>()?;
        let mut bases = 0;
        let mut sum = 0.0;
        for val in interval {
            let num_bases = val.end - val.start;
            bases += num_bases;
            sum += f64::from(num_bases) * f64::from(val.value);
        }
        let mean = if bases == 0 {
            0.0
        } else {
            sum / f64::from(bases)
        };
        bedoutwriter.write_fmt(format_args!("{}\t{}\t{}\t{}\t{:.3}\n", chrom, start, end, rest, mean))?;
    }
    Ok(())
}

fn main() -> io::Result<()> {
    let mut args = std::env::args();
    args.next();
    let bedinpath = args.next().expect("Must pass a bed input path.");
    let bigwigpath = args.next().expect("Must pass a bigwig input file.");
    let bedoutpath = args.next().expect("Must pass a bed output file.");
    println!("Args: {} {} {}", bedinpath, bigwigpath, bedoutpath);

    let inbed = File::open(bedinpath)?;
    let inbigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
    let outbed = File::create(bedoutpath)?;
    write(inbed, inbigwig, outbed)?;

    Ok(())
}
