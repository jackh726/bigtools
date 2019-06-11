#![feature(async_await)]

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};

use clap::{App, Arg};

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
        let matches = App::new("BigWigInfo")
        .arg(Arg::with_name("bedin")
                .help("The input bed file")
                .index(1)
                .required(true)
            )
        .arg(Arg::with_name("bigwig")
                .help("The input bigwig file")
                .index(2)
                .required(true)
            )
        .arg(Arg::with_name("output")
                .help("The output bed file")
                .index(3)
                .required(true)
            )
        .get_matches();

    let bedinpath = matches.value_of("bedin").unwrap().to_owned();
    let bigwigpath = matches.value_of("bigwig").unwrap().to_owned();
    let bedoutpath = matches.value_of("output").unwrap().to_owned();

    let inbed = File::open(bedinpath)?;
    let inbigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
    let outbed = File::create(bedoutpath)?;
    write(inbed, inbigwig, outbed)?;

    Ok(())
}
