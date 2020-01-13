use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};

use clap::{App, Arg, SubCommand};

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

    while let Some(line) = bedstream.read()? {
        let mut split = line.splitn(4, '\t');
        let chrom = split.next().expect("Missing chrom");
        let start = split.next().expect("Missing start").parse::<u32>().unwrap();
        let end = split.next().expect("Missing end").parse::<u32>().unwrap();
        let _ = split.next();
        let interval = b
            .get_interval(chrom, start, end);
        let interval = match interval {
            Ok(interval) => interval,
            Err(e) => {
                eprintln!("An error occured when intersecting {}:{}-{}: {}", chrom, start, end, e);
                continue;
            }
        };

        for raw_val in interval {
            let val = match raw_val {
                Ok(val) => val,
                Err(e) => {
                    eprintln!("An error occured when intersecting {}:{}-{}: {}", chrom, start, end, e);
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
            SubCommand::with_name("intersect")
                .about("Intersect all entries of a bed with a bigBed")
                .arg(
                    Arg::with_name("a")
                        .short("a")
                        .help("Each entry in this bed is compared against b for overlaps.")
                        .takes_value(true)
                        .required(true)
                )
                .arg(
                    Arg::with_name("b")
                        .short("b")
                        .help("Each entry in a will be compared against this bigBed for overlaps.")
                        .takes_value(true)
                        .required(true)
                )
        )
        .get_matches();

    match matches.subcommand() {
        ("intersect", Some(matches)) => {
            eprintln!("---BigTools intersect---");

            let apath = matches.value_of("a").unwrap().to_owned();
            let bpath = matches.value_of("b").unwrap().to_owned();

            let b = BigBedRead::from_file_and_attach(bpath)?;

            intersect(apath, b, IntersectOptions {})?;
        }
        (subcommand, _) => {
            panic!("BUG: unhandled subcommand: {}", subcommand);
        }
    }
    Ok(())
}
