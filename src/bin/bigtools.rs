use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};

use clap::{App, Arg};

use bigtools::bigwig::{BigBedRead, BigBedReadAttachError};
use bigtools::seekableread::{Reopen, SeekableRead};
use bigtools::streaming_linereader::StreamingLineReader;

struct IntersectOptions {}

fn intersect<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    apath: String,
    mut b: BigBedRead<R, S>,
    _options: IntersectOptions,
) -> io::Result<()> {
    let bedin = File::open(apath)?;
    let mut bedstream = StreamingLineReader::new(BufReader::with_capacity(64 * 1024, bedin));

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut bedoutwriter = BufWriter::with_capacity(64 * 1024, handle);

    while let Some(line) = bedstream.read() {
        let line = line?;
        let mut split = line.trim_end().splitn(4, '\t');
        let chrom = split.next().ok_or(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Missing chrom: {}", line),
        ))?;
        let start = split
            .next()
            .ok_or(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Missing start: {}", line),
            ))?
            .parse::<u32>()
            .map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid start: {:}", line),
                )
            })?;
        let end = split
            .next()
            .ok_or(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Missing end: {}", line),
            ))?
            .parse::<u32>()
            .map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid end: {:}", line),
                )
            })?;
        let interval = b.get_interval(chrom, start, end);
        let interval = match interval {
            Ok(interval) => interval,
            Err(e) => {
                eprintln!(
                    "An error occured when intersecting {}:{}-{}: {}",
                    chrom, start, end, e
                );
                continue;
            }
        };

        for raw_val in interval {
            let val = match raw_val {
                Ok(val) => val,
                Err(e) => {
                    eprintln!(
                        "An error occured when intersecting {}:{}-{}: {}",
                        chrom, start, end, e
                    );
                    continue;
                }
            };
            bedoutwriter.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\n",
                chrom, val.start, val.end, val.rest
            ))?;
        }
    }

    Ok(())
}

fn main() -> Result<(), BigBedReadAttachError> {
    let matches = App::new("BigTools")
        .subcommand(
            App::new("intersect")
                .about("Intersect all entries of a bed with a bigBed")
                .arg(
                    Arg::new("a")
                        .short('a')
                        .about("Each entry in this bed is compared against b for overlaps.")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("b")
                        .short('b')
                        .about("Each entry in a will be compared against this bigBed for overlaps.")
                        .takes_value(true)
                        .required(true),
                ),
        )
        .get_matches();

    match matches.subcommand() {
        Some(("intersect", matches)) => {
            eprintln!("---BigTools intersect---");

            let apath = matches.value_of("a").unwrap().to_owned();
            let bpath = matches.value_of("b").unwrap().to_owned();

            let b = BigBedRead::from_file_and_attach(bpath)?;

            intersect(apath, b, IntersectOptions {})?;
        }
        None => {
            eprintln!("No command. Use bigtools -.about to see.about.");
        }
        Some((subcommand, _)) => {
            panic!("BUG: unhandled subcommand: {}", subcommand);
        }
    }
    Ok(())
}
