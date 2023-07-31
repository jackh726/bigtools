use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};

use bigtools::{BigWigRead, BigWigReadAttachError};
use clap::{Arg, Command};

use bigtools::bbi::{BigBedRead, BigBedReadAttachError};
use bigtools::utils::reopen::SeekableRead;
use bigtools::utils::streaming_linereader::StreamingLineReader;

struct IntersectOptions {}

fn intersect<R: SeekableRead + 'static>(
    apath: String,
    mut b: BigBedRead<R>,
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
        let chrom = split.next().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Missing chrom: {}", line),
            )
        })?;
        let start = split
            .next()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Missing start: {}", line),
                )
            })?
            .parse::<u32>()
            .map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid start: {:}", line),
                )
            })?;
        let end = split
            .next()
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, format!("Missing end: {}", line))
            })?
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

fn chromintersect(apath: String, bpath: String, outpath: String) -> io::Result<()> {
    let chroms = match BigWigRead::open_file(&bpath) {
        Ok(bigwig) => bigwig.info.chrom_info,
        Err(BigWigReadAttachError::NotABigWig) => match BigBedRead::open_file(bpath) {
            Ok(bigbed) => bigbed.info.chrom_info,
            Err(BigBedReadAttachError::NotABigBed) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Only bigWigs and bigBeds are supported as `b` files."),
                ));
            }
            Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, format!("{}", e))),
        },
        Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, format!("{}", e))),
    };
    let chroms = HashSet::from_iter(chroms.into_iter().map(|c| c.name));

    fn write<T: Write>(
        chroms: HashSet<String>,
        apath: String,
        mut bedoutwriter: BufWriter<T>,
    ) -> io::Result<()> {
        let bedin = File::open(apath)?;
        let mut bedstream = StreamingLineReader::new(BufReader::with_capacity(64 * 1024, bedin));

        while let Some(line) = bedstream.read() {
            let line = line?;
            let mut split = line.trim_end().splitn(4, '\t');
            let chrom = split.next().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Missing chrom: {}", line),
                )
            })?;
            if chroms.contains(chrom) {
                bedoutwriter.write(line.as_bytes())?;
            }
        }
        Ok(())
    }

    if outpath == "-" {
        let stdout = io::stdout();
        let handle = stdout.lock();
        let bedoutwriter = BufWriter::with_capacity(64 * 1024, handle);
        write(chroms, apath, bedoutwriter)?;
    } else {
        let bedout = File::create(outpath)?;
        let bedoutwriter = BufWriter::with_capacity(64 * 1024, bedout);
        write(chroms, apath, bedoutwriter)?;
    }

    Ok(())
}

fn main() -> Result<(), BigBedReadAttachError> {
    let matches = Command::new("BigTools")
        .subcommand(
            Command::new("intersect")
                .about("Intersect all entries of a bed with a bigBed")
                .arg(
                    Arg::new("a")
                        .short('a')
                        .help("Each entry in this bed is compared against b for overlaps.")
                        .num_args(1)
                        .required(true),
                )
                .arg(
                    Arg::new("b")
                        .short('b')
                        .help("Each entry in a will be compared against this bigBed for overlaps.")
                        .num_args(1)
                        .required(true),
                ),
        )
        .subcommand(
            Command::new("chromintersect")
                .about("Create a new file of the same type, containing only data from `a` with chromosomes from `b`")
                .arg(
                    Arg::new("a")
                        .short('a')
                        .help("The file to take data from (currently supports: bed)")
                        .num_args(1)
                        .required(true),
                )
                .arg(
                    Arg::new("b")
                        .short('b')
                        .help("The file to take reference chromosomes from (currently supports: bigWig or bigBed)")
                        .num_args(1)
                        .required(true),
                )
                .arg(
                    Arg::new("out")
                        .short('o')
                        .help("The name of the output file (or - for stdout). Outputted in same format as `a`")
                        .num_args(1)
                        .required(true),
                ),
        )
        .get_matches();

    match matches.subcommand() {
        Some(("intersect", matches)) => {
            eprintln!("---BigTools intersect---");

            let apath = matches.get_one::<String>("a").unwrap().to_owned();
            let bpath = matches.get_one::<String>("b").unwrap().to_owned();

            let b = BigBedRead::open_file(bpath)?;

            intersect(apath, b, IntersectOptions {})?;
        }
        Some(("chromintersect", matches)) => {
            eprintln!("---BigTools chromintersect---");

            let apath = matches.get_one::<String>("a").unwrap().to_owned();
            let bpath = matches.get_one::<String>("b").unwrap().to_owned();
            let outpath = matches.get_one::<String>("out").unwrap().to_owned();

            chromintersect(apath, bpath, outpath)?;
        }
        None => {
            eprintln!("No command. Use bigtools -help to see help.");
        }
        Some((subcommand, _)) => {
            panic!("BUG: unhandled subcommand: {}", subcommand);
        }
    }
    Ok(())
}
