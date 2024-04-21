//! Utilities for parsing a bed file.
//!
//! There are roughly two layers of abstraction here, each with a different purpose.
//!
//! The first layer of abstraction is enscapsulated in the `StreamingBedValues` trait. Briefly,
//! implementors of this trait return "raw" bed-like data. This is the chromosome (as a `&str`) and
//! data specific to each type of bed
//!
//! The second layer of abstraction (`BedParser`) manages the state information for when values switch
//! from one chromosome to another. The is important because bigwig/bigbed writing is "chunked" by chromosome.

use std::io::{self, BufRead, BufReader, Read};

use thiserror::Error;

use crate::bbi::{BedEntry, Value};
use crate::utils::streaming_linereader::StreamingLineReader;

pub fn parse_bed<'a>(s: &'a str) -> Option<Result<(&'a str, BedEntry), BedValueError>> {
    let mut split = s.splitn(4, '\t');
    let chrom = match split.next() {
        Some(chrom) => chrom,
        None => return None,
    };
    let res = (|| {
        let s = split
            .next()
            .ok_or_else(|| BedValueError::InvalidInput(format!("Missing start: {:}", s)))?;
        let start = s
            .parse::<u32>()
            .map_err(|_| BedValueError::InvalidInput(format!("Invalid start: {:}", s)))?;
        let s = split
            .next()
            .ok_or_else(|| BedValueError::InvalidInput(format!("Missing end: {:}", s)))?;
        let end = s
            .parse::<u32>()
            .map_err(|_| BedValueError::InvalidInput(format!("Invalid end: {:}", s)))?;
        let rest = split.next().unwrap_or("").to_string();
        Ok((start, end, rest))
    })();
    match res {
        Err(e) => Some(Err(e)),
        Ok((start, end, rest)) => Some(Ok((chrom, BedEntry { start, end, rest }))),
    }
}

pub fn parse_bedgraph<'a>(s: &'a str) -> Option<Result<(&'a str, Value), BedValueError>> {
    let mut split = s.splitn(5, '\t');
    let chrom = match split.next() {
        Some(chrom) => chrom,
        None => return None,
    };
    let res = (|| {
        let s = split
            .next()
            .ok_or_else(|| BedValueError::InvalidInput(format!("Missing start: {:}", s)))?;
        let start = s
            .parse::<u32>()
            .map_err(|_| BedValueError::InvalidInput(format!("Invalid start: {:}", s)))?;
        let s = split
            .next()
            .ok_or_else(|| BedValueError::InvalidInput(format!("Missing end: {:}", s)))?;
        let end = s
            .parse::<u32>()
            .map_err(|_| BedValueError::InvalidInput(format!("Invalid end: {:}", s)))?;
        let s = split
            .next()
            .ok_or_else(|| BedValueError::InvalidInput(format!("Missing value: {:}", s)))?;
        let value = s
            .parse::<f32>()
            .map_err(|_| BedValueError::InvalidInput(format!("Invalid value: {:}", s)))?;
        Ok((start, end, value))
    })();
    match res {
        Err(e) => Some(Err(e)),
        Ok((start, end, value)) => Some(Ok((chrom, Value { start, end, value }))),
    }
}

// FIXME: can replace with this with just a simple `LendingIterator`
/// Essentially a combined lending iterator over the chrom (&str) and remaining
/// values of bed-like data
pub trait StreamingBedValues {
    type Value;

    fn next(&mut self) -> Option<Result<(&str, Self::Value), BedValueError>>;
}

#[derive(Error, Debug)]
pub enum BedValueError {
    #[error("{}", .0)]
    InvalidInput(String),
    #[error("Error occurred: {}", .0)]
    IoError(#[from] io::Error),
}

pub type Parser<V> = for<'a> fn(&'a str) -> Option<Result<(&'a str, V), BedValueError>>;

/// Parses a bed-like file
pub struct BedFileStream<V, B> {
    pub bed: StreamingLineReader<B>,
    pub parse: Parser<V>,
}

impl<R: Read> BedFileStream<BedEntry, BufReader<R>> {
    pub fn from_bed_file(file: R) -> BedFileStream<BedEntry, BufReader<R>> {
        BedFileStream {
            bed: StreamingLineReader::new(BufReader::new(file)),
            parse: parse_bed,
        }
    }
}

impl<R: Read> BedFileStream<Value, BufReader<R>> {
    pub fn from_bedgraph_file(file: R) -> BedFileStream<Value, BufReader<R>> {
        BedFileStream {
            bed: StreamingLineReader::new(BufReader::new(file)),
            parse: parse_bedgraph,
        }
    }
}

impl<V, B: BufRead> StreamingBedValues for BedFileStream<V, B> {
    type Value = V;

    fn next(&mut self) -> Option<Result<(&str, Self::Value), BedValueError>> {
        let line = match self.bed.read()? {
            Ok(line) => line.trim_end(),
            Err(e) => return Some(Err(e.into())),
        };
        match (self.parse)(line) {
            None => None,
            Some(Ok(v)) => Some(Ok(v)),
            Some(Err(e)) => Some(Err(e.into())),
        }
    }
}

// Wraps a bed-like Iterator
pub struct BedIteratorStream<V, I> {
    pub(crate) iter: I,
    pub(crate) curr: Option<(String, V)>,
}

impl<
        V: Clone,
        E: Into<BedValueError>,
        C: Into<String> + for<'a> PartialEq<&'a str>,
        I: Iterator<Item = Result<(C, V), E>>,
    > StreamingBedValues for BedIteratorStream<V, I>
{
    type Value = V;

    fn next(&mut self) -> Option<Result<(&str, V), BedValueError>> {
        use std::ops::Deref;
        self.curr = match (self.curr.take(), self.iter.next()?) {
            (_, Err(e)) => return Some(Err(e.into())),
            (Some(c), Ok(v)) => {
                if v.0 == &c.0 {
                    Some((c.0, v.1))
                } else {
                    Some((v.0.into(), v.1))
                }
            }
            (None, Ok(v)) => Some((v.0.into(), v.1)),
        };
        self.curr.as_ref().map(|v| Ok((v.0.deref(), v.1.clone())))
    }
}

// Wraps a bed-like Iterator
pub struct BedInfallibleIteratorStream<V, I> {
    pub(crate) iter: I,
    pub(crate) curr: Option<(String, V)>,
}

impl<V: Clone, C: Into<String> + for<'a> PartialEq<&'a str>, I: Iterator<Item = (C, V)>>
    StreamingBedValues for BedInfallibleIteratorStream<V, I>
{
    type Value = V;

    fn next(&mut self) -> Option<Result<(&str, V), BedValueError>> {
        use std::ops::Deref;
        self.curr = match (self.curr.take(), self.iter.next()?) {
            (Some(c), v) => {
                if v.0 == &c.0 {
                    Some((c.0, v.1))
                } else {
                    Some((v.0.into(), v.1))
                }
            }
            (None, v) => Some((v.0.into(), v.1)),
        };
        self.curr.as_ref().map(|v| Ok((v.0.deref(), v.1.clone())))
    }
}
