use std::collections::HashSet;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};

use bigtools::{BBIRead, GenericBBIRead};
use clap::{Parser, Subcommand};

use bigtools::utils::reopen::SeekableRead;
use bigtools::utils::streaming_linereader::StreamingLineReader;
use bigtools::BigBedRead;

#[derive(Debug, Subcommand)]
enum SubCommands {
    Intersect {
        /// Each entry in this bed is compared against `b` for overlaps.
        a: String,

        /// Each entry in `a` will be compared against this bigBed for overlaps.
        b: String,
    },
    ChromIntersect {
        /// The file to take data from (currently supports: bed)
        a: String,

        /// The file to take reference chromosomes from (currently supports: bigWig or bigBed)
        b: String,

        /// The name of the output file (or - for stdout). Outputted in same format as `a`.
        out: String,
    },
}

#[derive(Debug, Parser)]
#[command(about = "BigTools", long_about = None, multicall = true)]
enum CliCommands {
    Bigtools {
        #[command(subcommand)]
        command: SubCommands,
    },
    #[command(flatten)]
    SubCommands(SubCommands),
}

struct IntersectOptions {}

fn intersect<R: SeekableRead + 'static>(
    apath: String,
    mut b: BigBedRead<R>,
    _options: IntersectOptions,
) -> io::Result<()> {
    let bedin = File::open(&apath)?;
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
    let chroms = match GenericBBIRead::open_file(&bpath) {
        Ok(b) => b.chroms().to_vec(),
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

fn main() -> Result<(), Box<dyn Error>> {
    let args = env::args_os();
    let mut args: Vec<_> = args.collect();
    args[0] = args[0].to_ascii_lowercase();
    let cli = CliCommands::parse_from(args);
    let command = match cli {
        CliCommands::Bigtools { command } => command,
        CliCommands::SubCommands(command) => command,
    };
    match command {
        SubCommands::Intersect { a, b } => {
            let b = BigBedRead::open_file(&b)?;

            intersect(a, b, IntersectOptions {})?;
        }
        SubCommands::ChromIntersect { a, b, out } => {
            chromintersect(a, b, out)?;
        }
    }
    Ok(())
}
