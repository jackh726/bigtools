use std::io::{self, BufRead};

use crate::bbiread::BBIReadError;
use crate::bigwig::BigWigRead;
use crate::utils::file::seekableread::{Reopen, SeekableRead};
use crate::utils::file::streaming_linereader::StreamingLineReader;

#[derive(Copy, Clone)]
pub enum Name {
    Interval,
    None,
    Column(usize),
}

pub struct BigWigAverageOverBedEntry {
    pub name: String,
    pub size: u32,
    pub bases: u32,
    pub sum: f64,
    pub mean0: f64,
    pub mean: f64,
}

pub fn stats_for_bed_item<R: Reopen<S>, S: SeekableRead>(
    name: Name,
    line: &str,
    bigwig: &mut BigWigRead<R, S>,
) -> Result<BigWigAverageOverBedEntry, BBIReadError> {
    let cols = line.trim_end().split('\t').collect::<Vec<_>>();
    if cols.len() < 3 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid bed: A minimum of 3 columns must be specified (chrom, start, end).",
        )
        .into());
    }

    let chrom = cols[0];
    let start = cols[1]
        .parse::<u32>()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid start: not an integer"))?;
    let end = cols[2]
        .parse::<u32>()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid end: not an integer"))?;
    let interval = bigwig
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

    let name = match name {
        Name::Column(col) => {
            match cols.get(col) {
                Some(v) => v.to_string(),
                None => return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid name column option. Number of columns ({}) is less than the value specified ({}).", cols.len(), col+1)
                ).into()),
            }
        }
        Name::Interval => format!("{}:{}-{}", chrom, start, end),
        Name::None => cols.join("\t"),
    };

    Ok(BigWigAverageOverBedEntry {
        name,
        size,
        bases,
        sum,
        mean0,
        mean,
    })
}

pub fn bigwig_average_over_bed<R: Reopen<S> + 'static, S: SeekableRead + 'static>(
    bed: impl BufRead,
    mut bigwig: BigWigRead<R, S>,
    name: Name,
) -> Result<impl Iterator<Item = Result<BigWigAverageOverBedEntry, BBIReadError>>, BBIReadError> {
    let mut bedstream = StreamingLineReader::new(bed);

    let mut error: bool = false;
    let iter = std::iter::from_fn(
        move || -> Option<Result<BigWigAverageOverBedEntry, BBIReadError>> {
            if error {
                return None;
            }
            let line = match bedstream.read()? {
                Err(e) => {
                    error = true;
                    return Some(Err(e.into()));
                }
                Ok(line) => line,
            };
            match stats_for_bed_item(name, line, &mut bigwig) {
                Err(e) => {
                    error = true;
                    Some(Err(e))
                }
                Ok(v) => Some(Ok(v)),
            }
        },
    );

    Ok(iter)
}
