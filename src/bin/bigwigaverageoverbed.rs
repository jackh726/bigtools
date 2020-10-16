use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use clap::{App, Arg};

use bigtools::bigwig::{BigWigRead, BigWigReadAttachError};
use bigtools::utils::misc::{bigwig_average_over_bed, BigWigAverageOverBedOptions, Name};

fn main() -> Result<(), BigWigReadAttachError> {
    let matches = App::new("BigWigAverageOverBed")
        .arg(Arg::new("bigwig")
                .about("The input bigwig file")
                .index(1)
                .required(true)
            )
        .arg(Arg::new("bedin")
                .about("The input bed file")
                .index(2)
                .required(true)
            )
        .arg(Arg::new("output")
                .about("The output bed file")
                .index(3)
                .required(true)
            )
        .arg(Arg::new("allcols")
                .short('a')
                .about("If set, the output will be a bed with all columns of the input, with additional states at the end. Otherwise, the ouput will be a tsv with the name (or the chrom, start, and end in the form `chrom:start-end` if names are not unique) followed by stats.")
            )
        .get_matches();

    let bigwigpath = matches.value_of("bigwig").unwrap();
    let bedinpath = matches.value_of("bedin").unwrap();
    let bedoutpath = matches.value_of("output").unwrap();

    let allcols = matches.is_present("allcols");

    let inbigwig = BigWigRead::from_file_and_attach(bigwigpath)?;
    let outbed = File::create(bedoutpath)?;
    let mut bedoutwriter = BufWriter::new(outbed);

    let name: Name = {
        if allcols {
            Name::None
        } else {
            const NAME_COL: usize = 4;
            const TEST_LINES: usize = 10;
            let reader = BufReader::new(File::open(bedinpath.clone())?);
            let mut lines = reader
                .lines()
                .take(TEST_LINES)
                .map(|line| -> io::Result<Option<String>> {
                    let l = line?;
                    let cols = l
                        .trim()
                        .split('\t')
                        .take((NAME_COL + 1).max(3))
                        .collect::<Vec<_>>();
                    if cols.len() < 3 {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Not enough columns. Expected >= 3.",
                        ));
                    }
                    Ok(cols.get(NAME_COL - 1).map(|s| s.to_string()))
                })
                .collect::<io::Result<Vec<_>>>()?;
            lines.sort();
            lines.dedup();
            if lines.len() != TEST_LINES {
                eprintln!(
                    "Name column ({}) is not unique. Using interval as the name.",
                    NAME_COL
                );
                Name::None
            } else {
                Name::Column(std::num::NonZeroUsize::new(NAME_COL).unwrap())
            }
        }
    };

    let options = BigWigAverageOverBedOptions { name, allcols };
    let bedin = BufReader::new(File::open(bedinpath)?);
    let iter = bigwig_average_over_bed(bedin, inbigwig, options)?;
    for entry in iter {
        let entry = entry?;
        let name = entry.name;
        let stats = format!(
            "{}\t{}\t{:.3}\t{:.3}\t{:.3}",
            entry.size, entry.bases, entry.sum, entry.mean0, entry.mean
        );
        write!(&mut bedoutwriter, "{}\t{}\n", name, stats)?
    }

    Ok(())
}
