use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};

use bigtools::utils::streaming_linereader::StreamingLineReader;
use clap::{App, Arg};

use bigtools::bigwig::BigWigRead;
use bigtools::utils::misc::{stats_for_bed_item, Name};

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("BigWigAverageOverBed")
        .arg(Arg::new("bigwig")
                .help("The input bigwig file")
                .index(1)
                .required(true)
            )
        .arg(Arg::new("bedin")
                .help("The input bed file")
                .index(2)
                .required(true)
            )
        .arg(Arg::new("output")
                .help("The output bed file")
                .index(3)
                .required(true)
            )
        .arg(Arg::new("namecol")
                .short('n')
                .help("Supports three types of options: `interval`, `none`, or a column number (one indexed). If `interval`, the name column in the output will be the interval in the form of `chrom:start-end`. If `none`, then all columns will be included in the output file. Otherwise, the one-indexed column will be used as the name. By default, column 4 is used as a name column.")
                .default_value("4")
            )
        .get_matches();

    let bigwigpath = matches.value_of("bigwig").unwrap();
    let bedinpath = matches.value_of("bedin").unwrap();
    let bedoutpath = matches.value_of("output").unwrap();

    let mut inbigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
    let outbed = File::create(bedoutpath)?;
    let mut bedoutwriter = BufWriter::new(outbed);

    let bedin = BufReader::new(File::open(bedinpath)?);
    let mut bedstream = StreamingLineReader::new(bedin);

    let name = match matches.value_of("namecol") {
        Some("interval") => Name::Interval,
        Some("none") => Name::None,
        Some(col) => {
            let col = col.parse::<usize>();
            match col {
                Ok(col) if col > 0 => Name::Column(col - 1),
                Ok(_) => return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid name column option. Column values are one-indexed, so should not be zero.",
                ).into()),
                Err(_) => return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid name column option. Allowed options are `interval`, `none`, or an one-indexed integer value for a given column.",
                ).into()),
            }
        }
        None => Name::Column(3),
    };

    loop {
        let line = match bedstream.read() {
            None => break,
            Some(Err(e)) => {
                return Err(e.into());
            }
            Some(Ok(line)) => line,
        };

        let entry = match stats_for_bed_item(name, line, &mut inbigwig) {
            Ok(stats) => stats,
            Err(e) => return Err(e.into()),
        };

        let stats = format!(
            "{}\t{}\t{:.3}\t{:.3}\t{:.3}",
            entry.size, entry.bases, entry.sum, entry.mean0, entry.mean
        );
        writeln!(&mut bedoutwriter, "{}\t{}", entry.name, stats)?
    }

    Ok(())
}
