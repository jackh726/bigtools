use std::io::BufRead;

use crate::bbi::BigWigRead;
use crate::bbiread::BBIReadError;
use crate::bed::bedparser::{parse_bed, BedValueError};
use crate::utils::file::streaming_linereader::StreamingLineReader;
use crate::{BBIFileRead, BedEntry};

use thiserror::Error;

#[derive(Copy, Clone, Debug)]
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

#[derive(Error, Debug)]
pub enum StatsError {
    #[error("{}", .0)]
    BBIReadError(#[from] BBIReadError),
    #[error("{}", .0)]
    InvalidNameCol(String),
}

pub fn stats_for_bed_item<R: BBIFileRead>(
    name: Name,
    chrom: &str,
    entry: BedEntry,
    bigwig: &mut BigWigRead<R>,
) -> Result<BigWigAverageOverBedEntry, StatsError> {
    let start = entry.start;
    let end = entry.end;

    let interval = bigwig
        .get_interval(chrom, start, end)?
        .collect::<Result<Vec<_>, _>>();
    let interval = match interval {
        Ok(i) => i,
        Err(e) => return Err(e.into()),
    };

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
        Name::Column(col) => match col {
            0 => chrom.to_string(),
            1 => start.to_string(),
            2 => end.to_string(),
            _ => {
                let mut cols = entry.rest.split('\t');
                let v = cols.nth(col - 3);
                match v {
                        Some(v) => v.to_string(),
                        None => return Err(StatsError::InvalidNameCol(format!("Invalid name column option. Number of columns ({}) is less than the value specified ({}).", entry.rest.split('\t').collect::<Vec<_>>().len() + 3, col+1))),
                    }
            }
        },
        Name::Interval => format!("{}:{}-{}", chrom, start, end),
        Name::None => format!("{}\t{}\t{}\t{}", chrom, start, end, entry.rest),
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

#[derive(Error, Debug)]
pub enum BigWigAverageOverBedError {
    #[error("{}", .0)]
    StatsError(#[from] StatsError),
    #[error("{}", .0)]
    BedValueError(#[from] BedValueError),
}

pub fn bigwig_average_over_bed<R: BBIFileRead>(
    bed: impl BufRead,
    mut bigwig: BigWigRead<R>,
    name: Name,
) -> impl Iterator<Item = Result<BigWigAverageOverBedEntry, BigWigAverageOverBedError>> {
    let mut bedstream = StreamingLineReader::new(bed);

    let mut error: bool = false;
    let iter = std::iter::from_fn(
        move || -> Option<Result<BigWigAverageOverBedEntry, BigWigAverageOverBedError>> {
            if error {
                return None;
            }
            let line = match bedstream.read()? {
                Err(e) => {
                    error = true;
                    return Some(Err(BedValueError::IoError(e).into()));
                }
                Ok(line) => line,
            };
            let (chrom, entry) = match parse_bed(line) {
                None => return None,
                Some(Err(e)) => {
                    error = true;
                    return Some(Err(e.into()));
                }
                Some(Ok(v)) => v,
            };

            match stats_for_bed_item(name, chrom, entry, &mut bigwig) {
                Err(e) => {
                    error = true;
                    Some(Err(e.into()))
                }
                Ok(v) => Some(Ok(v)),
            }
        },
    );

    iter
}
