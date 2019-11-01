use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use clap::{App, Arg};

use bigtools::bigwig::{BigWigRead, BigWigReadAttachError};
use bigtools::seekableread::{Reopen, SeekableRead};
use bigtools::streaming_linereader::StreamingLineReader;

struct Options {
    simple: bool,
}

fn write<R: Reopen<S> + 'static, S: SeekableRead + 'static>(bedinpath: String, mut bigwigin: BigWigRead<R, S>, bedout: File, options: Options) -> io::Result<()> {
    let uniquenames = {
        if !options.simple {
            true
        } else {
            let reader = BufReader::new(File::open(bedinpath.clone())?); 
            let mut lines = reader
                .lines()
                .take(10)
                .map(|line| -> io::Result<Option<String>>{
                    let l = line?;
                    let mut split = l.splitn(5, '\t');
                    let _chrom = split.next().ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData));
                    let _start = split.next().ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData));
                    let _end = split.next().ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData));
                    let name = split.next();
                    Ok(match name {
                        None => None,
                        Some(name) => Some(name.to_string()),
                    })
                })
                .collect::<io::Result<Vec<_>>>()?;
            lines.sort();
            lines.dedup();
            lines.len() == 10
        }
    };

    let bedin = File::open(bedinpath)?;
    let mut bedstream = StreamingLineReader::new(BufReader::new(bedin));
    let mut bedoutwriter = BufWriter::new(bedout);

    while let Some(line) = bedstream.read()? {
        let mut split = line.splitn(5, '\t');
        let chrom = split.next().expect("Missing chrom");
        let start = split.next().expect("Missing start").parse::<u32>().unwrap();
        let end = split.next().expect("Missing end").parse::<u32>().unwrap();
        let name = split.next();
        let rest = split.next();
        let interval = bigwigin.get_interval(chrom, start, end)?.collect::<Result<Vec<_>, _>>()?;
        let mut bases = 0;
        let mut sum = 0.0;
        for val in interval {
            let num_bases = val.end - val.start;
            bases += num_bases;
            sum += f64::from(num_bases) * f64::from(val.value);
        }
        let size = end - start;
        let mean0 = sum / f64::from(size);
        let mean = if bases == 0 {
            0.0
        } else {
            sum / f64::from(bases)
        };
        let stats = format!("{}\t{}\t{:.3}\t{:.3}\t{:.3}", size, bases, sum, mean0, mean);
        if options.simple {
            if uniquenames {
                bedoutwriter.write_fmt(format_args!("{}\t{}\n", name.expect("Bad bed format (no name)."), stats))?;
            } else {
                bedoutwriter.write_fmt(format_args!("{}\t{}\t{}\t{}\n", chrom, start, end, stats))?;
            }
        } else {
            let last = match name {
                Some(name) => match rest {
                    Some(rest) => format!("{}\t{}", name, rest.trim_end()),
                    None => name.trim_end().to_string(),
                },
                None => String::from(""),
            };
            bedoutwriter.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\n", chrom, start, end, last, stats))?;
        }
    }
    Ok(())
}

fn main() -> Result<(), BigWigReadAttachError> {
        let matches = App::new("BigWigAverageOverBed")
        .arg(Arg::with_name("bigwig")
                .help("The input bigwig file")
                .index(1)
                .required(true)
            )
        .arg(Arg::with_name("bedin")
                .help("The input bed file")
                .index(2)
                .required(true)
            )
        .arg(Arg::with_name("output")
                .help("The output bed file")
                .index(3)
                .required(true)
            )
        .arg(Arg::with_name("simple")
                .short("s")
                .help("If set, the output bed will either print the name followed by stats (if the fourth column in the first 10 entries), or the chrom,start,end followed by stats. Any other columns in the input bed will be ignored.")
            )
        .get_matches();

    let bigwigpath = matches.value_of("bigwig").unwrap().to_owned();
    let bedinpath = matches.value_of("bedin").unwrap().to_owned();
    let bedoutpath = matches.value_of("output").unwrap().to_owned();

    let simple = matches.is_present("simple");

    let inbigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
    let outbed = File::create(bedoutpath)?;
    let options = Options {
        simple,
    };
    write(bedinpath, inbigwig, outbed, options)?;

    Ok(())
}
