use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};

use bigtools::utils::cli::bedgraphtobigwig::{bedgraphtobigwig, BedGraphToBigWigArgs};
use bigtools::utils::cli::bedtobigbed::{bedtobigbed, BedToBigBedArgs};
use bigtools::utils::cli::bigbedinfo::{bigbedinfo, BigBedInfoArgs};
use bigtools::utils::cli::bigbedtobed::{bigbedtobed, BigBedToBedArgs};
use bigtools::utils::cli::bigwigaverageoverbed::{bigwigaverageoverbed, BigWigAverageOverBedArgs};
use bigtools::utils::cli::bigwiginfo::{bigwiginfo, BigWigInfoArgs};
use bigtools::utils::cli::bigwigmerge::{bigwigmerge, BigWigMergeArgs};
use bigtools::utils::cli::bigwigtobedgraph::{bigwigtobedgraph, BigWigToBedGraphArgs};
use bigtools::utils::cli::bigwigvaluesoverbed::{bigwigvaluesoverbed, BigWigValuesOverBedArgs};
use bigtools::utils::cli::compat_args;
use bigtools::{BBIRead, GenericBBIRead};
use clap::{Args, Parser, Subcommand};

use bigtools::utils::reopen::SeekableRead;
use bigtools::utils::streaming_linereader::StreamingLineReader;
use bigtools::BigBedRead;

#[derive(Clone, Debug, PartialEq, Args)]
struct IntersectArgs {
    /// Each entry in this bed is compared against `b` for overlaps.
    a: String,

    /// Each entry in `a` will be compared against this bigBed for overlaps.
    b: String,
}

#[derive(Clone, Debug, PartialEq, Args)]
struct ChromIntersectArgs {
    /// The file to take data from (currently supports: bed)
    a: String,

    /// The file to take reference chromosomes from (currently supports: bigWig or bigBed)
    b: String,

    /// The name of the output file (or - for stdout). Outputted in same format as `a`.
    out: String,
}

#[derive(Clone, Debug, PartialEq, Subcommand)]
#[command(version)]
enum SubCommands {
    Intersect {
        #[command(flatten)]
        args: IntersectArgs,
    },
    #[command(name = "chromintersect", version)]
    ChromIntersect {
        #[command(flatten)]
        args: ChromIntersectArgs,
    },
    #[command(name = "bedgraphtobigwig", version)]
    BedGraphToBigWig {
        #[command(flatten)]
        args: BedGraphToBigWigArgs,
    },
    #[command(name = "bedtobigbed", version)]
    BedToBigBed {
        #[command(flatten)]
        args: BedToBigBedArgs,
    },
    #[command(name = "bigbedinfo", version)]
    BigBedInfo {
        #[command(flatten)]
        args: BigBedInfoArgs,
    },
    #[command(name = "bigbedtobed", version)]
    BigBedToBed {
        #[command(flatten)]
        args: BigBedToBedArgs,
    },
    #[command(name = "bigwigaverageoverbed", version)]
    BigWigAverageOverBed {
        #[command(flatten)]
        args: BigWigAverageOverBedArgs,
    },
    #[command(name = "bigwiginfo", version)]
    BigWigInfo {
        #[command(flatten)]
        args: BigWigInfoArgs,
    },
    #[command(name = "bigwigmerge", version)]
    BigWigMerge {
        #[command(flatten)]
        args: BigWigMergeArgs,
    },
    #[command(name = "bigwigtobedgraph", version)]
    BigWigToBedGraph {
        #[command(flatten)]
        args: BigWigToBedGraphArgs,
    },
    #[command(name = "bigwigvaluesoverbed", version)]
    BigWigValuesOverBed {
        #[command(flatten)]
        args: BigWigValuesOverBedArgs,
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
) -> Result<(), Box<dyn Error>> {
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

fn chromintersect(apath: String, bpath: String, outpath: String) -> Result<(), Box<dyn Error>> {
    let chroms = match GenericBBIRead::open_file(&bpath) {
        Ok(b) => b.chroms().to_vec(),
        Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, format!("{}", e)).into()),
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
    let cli = CliCommands::parse_from(compat_args(std::env::args_os()));
    let command = match cli {
        CliCommands::Bigtools { command } => command,
        CliCommands::SubCommands(command) => command,
    };
    match command {
        SubCommands::Intersect {
            args: IntersectArgs { a, b },
        } => {
            let b = BigBedRead::open_file(&b)?;

            intersect(a, b, IntersectOptions {})
        }
        SubCommands::ChromIntersect {
            args: ChromIntersectArgs { a, b, out },
        } => chromintersect(a, b, out),
        SubCommands::BedGraphToBigWig { args } => bedgraphtobigwig(args),
        SubCommands::BedToBigBed { args } => bedtobigbed(args),
        SubCommands::BigBedInfo { args } => bigbedinfo(args),
        SubCommands::BigBedToBed { args } => bigbedtobed(args),
        SubCommands::BigWigAverageOverBed { args } => {
            match bigwigaverageoverbed(args) {
                Ok(_) => {}
                Err(e) => {
                    // Returns `dyn Error + Send + Sync`, which can't be converted to `dyn Error`
                    // Just print and exit
                    eprintln!("{e}");
                    std::process::exit(1);
                }
            }

            Ok(())
        }
        SubCommands::BigWigInfo { args } => bigwiginfo(args),
        SubCommands::BigWigMerge { args } => bigwigmerge(args),
        SubCommands::BigWigToBedGraph { args } => bigwigtobedgraph(args),
        SubCommands::BigWigValuesOverBed { args } => bigwigvaluesoverbed(args),
    }
}
