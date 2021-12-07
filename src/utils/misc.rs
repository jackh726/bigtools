use std::io::{self, BufRead};
use std::num::NonZeroUsize;

use crate::bigwig::BigWigRead;
use crate::seekableread::{Reopen, SeekableRead};
use crate::streaming_linereader::StreamingLineReader;

pub enum Name {
    None,
    Column(NonZeroUsize),
}

pub struct BigWigAverageOverBedOptions {
    pub name: Name,
    pub allcols: bool,
}

pub struct BigWigAverageOverBedEntry {
    pub name: String,
    pub size: u32,
    pub bases: u32,
    pub sum: f64,
    pub mean0: f64,
    pub mean: f64,
}

pub fn bigwig_average_over_bed<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    bed: impl BufRead,
    mut bigwigin: BigWigRead<R, S>,
    options: BigWigAverageOverBedOptions,
) -> io::Result<impl Iterator<Item = io::Result<BigWigAverageOverBedEntry>>> {
    let mut bedstream = StreamingLineReader::new(bed);

    let mut error: bool = false;
    let iter = std::iter::from_fn(move || -> Option<io::Result<BigWigAverageOverBedEntry>> {
        if error {
            return None;
        }
        let line = match bedstream.read()? {
            Err(e) => {
                error = true;
                return Some(Err(e));
            }
            Ok(line) => line,
        };
        match (|| {
            let take_col = match options.name {
                Name::None => 0,
                Name::Column(col) => col.get(),
            };
            // We split by one more than we will be reading because we don't want "extra" data in the columns we check
            let cols = line
                .trim_end()
                .splitn(3.max(take_col) + 1, '\t')
                .collect::<Vec<_>>();
            if cols.len() < 3 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid bed: A minimum of 3 columns must be specified (chrom, start, end).",
                ));
            }
            if take_col != 0 && cols.len() < take_col {
                return Err(io::Error::new(io::ErrorKind::InvalidData, format!("Invalid bed: The column used for the name ({}) extends past the number of columns.", take_col)));
            }
            let interval_cols = &cols[0..=3.max(take_col)];
            let chrom = interval_cols[0];
            let start = interval_cols[1].parse::<u32>().map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "Invalid start: not an integer")
            })?;
            let end = interval_cols[2].parse::<u32>().map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "Invalid end: not an integer")
            })?;
            let interval = bigwigin
                .get_interval(chrom, start, end)?
                .collect::<Result<Vec<_>, _>>()?;

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

            let name = if options.allcols {
                cols.join("\t")
            } else {
                match options.name {
                    Name::None => format!("{}:{}-{}", chrom, start, end),
                    Name::Column(col) => interval_cols[col.get() - 1].to_string(),
                }
            };

            Ok(BigWigAverageOverBedEntry {
                name,
                size,
                bases,
                sum,
                mean0,
                mean,
            })
        })() {
            Err(e) => {
                error = true;
                Some(Err(e))
            }
            Ok(v) => Some(Ok(v)),
        }
    });

    Ok(iter)
}
