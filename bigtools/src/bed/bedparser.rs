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
use std::sync::Arc;

use crossbeam_utils::atomic::AtomicCell;
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
    iter: I,
    curr: Option<(String, V)>,
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
    iter: I,
    curr: Option<(String, V)>,
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

/// A wrapper for "bed-like" data
pub struct BedParser<S: StreamingBedValues> {
    state: Arc<AtomicCell<Option<BedParserState<S>>>>,
}

/// Defines the internal states of bed parsing
pub(crate) enum StateValue<V> {
    // No value has been loaded yet
    Empty,
    // A value has been loaded without error
    // Contains the current chromosome and the value.
    Value(String, V),
    // A previously loaded value was taken.
    // Contains the current chromosome.
    EmptyValue(String),
    // A new chromsome has been loaded
    DiffChrom(String, V),
    // An error has been seen
    Error(BedValueError),
    // We are done, either because we have run out of values or because of an error
    Done,
}

impl<V> std::fmt::Debug for StateValue<V> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Empty => write!(f, "Empty"),
            Self::Value(arg0, _) => f.debug_tuple("Value").field(arg0).finish(),
            Self::EmptyValue(arg0) => f.debug_tuple("EmptyValue").field(arg0).finish(),
            Self::DiffChrom(arg0, _) => f.debug_tuple("DiffChrom").field(arg0).finish(),
            Self::Error(arg0) => f.debug_tuple("Error").field(arg0).finish(),
            Self::Done => write!(f, "Done"),
        }
    }
}

impl<V> StateValue<V> {
    fn take_error(&mut self) -> Option<BedValueError> {
        let v = std::mem::replace(self, StateValue::Done);
        let ret;
        (*self, ret) = match v {
            StateValue::Error(e) => (StateValue::Done, Some(e)),
            s => (s, None),
        };
        ret
    }

    fn active_chrom(&self) -> Option<&String> {
        match self {
            StateValue::Empty => None,
            StateValue::Value(c, _) => Some(c),
            StateValue::EmptyValue(c) => Some(c),
            StateValue::DiffChrom(c, _) => Some(c),
            StateValue::Error(_) => None,
            StateValue::Done => None,
        }
    }
}

#[derive(Debug)]
pub(crate) struct BedParserState<S: StreamingBedValues> {
    stream: S,
    pub(crate) state_value: StateValue<S::Value>,
}

impl<S: StreamingBedValues> BedParser<S> {
    pub fn new(stream: S) -> Self {
        let state = BedParserState {
            stream,
            state_value: StateValue::Empty,
        };
        BedParser {
            state: Arc::new(AtomicCell::new(Some(state))),
        }
    }
}

impl<R: Read> BedParser<BedFileStream<BedEntry, BufReader<R>>> {
    pub fn from_bed_file(file: R) -> Self {
        BedParser::new(BedFileStream {
            bed: StreamingLineReader::new(BufReader::new(file)),
            parse: parse_bed,
        })
    }
}

impl<R: Read> BedParser<BedFileStream<Value, BufReader<R>>> {
    pub fn from_bedgraph_file(file: R) -> Self {
        BedParser::new(BedFileStream {
            bed: StreamingLineReader::new(BufReader::new(file)),
            parse: parse_bedgraph,
        })
    }
}

impl<
        V: Clone,
        E: Into<BedValueError>,
        C: Into<String> + for<'a> PartialEq<&'a str>,
        I: Iterator<Item = Result<(C, V), E>>,
    > BedParser<BedIteratorStream<V, I>>
{
    pub fn wrap_iter(iter: I) -> Self {
        BedParser::new(BedIteratorStream { iter, curr: None })
    }
}

impl<V: Clone, C: Into<String> + for<'a> PartialEq<&'a str>, I: Iterator<Item = (C, V)>>
    BedParser<BedInfallibleIteratorStream<V, I>>
{
    pub fn wrap_infallible_iter(iter: I) -> Self {
        BedParser::new(BedInfallibleIteratorStream { iter, curr: None })
    }
}

impl<S: StreamingBedValues> BedParser<S> {
    // This is *valid* to call multiple times for the same chromosome (assuming the
    // `BedChromData` has been dropped), since calling this function doesn't
    // actually advance the state (it will only set `next_val` if it currently is none).
    pub fn next_chrom(&mut self) -> Option<Result<(String, BedChromData<S>), BedValueError>> {
        let mut state = self.state.swap(None).expect("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
        state.load_state(true);
        let error = state.state_value.take_error();
        let chrom = state.state_value.active_chrom().cloned();
        self.state.swap(Some(state));

        if let Some(e) = error {
            return Some(Err(e));
        }

        match chrom {
            Some(chrom) => {
                let group = BedChromData {
                    state: self.state.clone(),
                    curr_state: None,
                    done: false,
                };
                Some(Ok((chrom.to_owned(), group)))
            }
            None => None,
        }
    }
}

impl<S: StreamingBedValues> BedParserState<S> {
    pub(crate) fn load_state(&mut self, switch_chrom: bool) {
        let state_value = std::mem::replace(&mut self.state_value, StateValue::Empty);
        self.state_value = match state_value {
            StateValue::Empty => match self.stream.next() {
                None => StateValue::Done,
                Some(Ok((chrom, val))) => StateValue::Value(chrom.to_owned(), val),
                Some(Err(err)) => StateValue::Error(err),
            },
            StateValue::Value(c, v) => StateValue::Value(c, v),
            StateValue::EmptyValue(prev_chrom) => match self.stream.next() {
                None => StateValue::Done,
                Some(Ok((chrom, val))) if switch_chrom || prev_chrom == chrom => {
                    StateValue::Value(prev_chrom, val)
                }
                Some(Ok((chrom, val))) => StateValue::DiffChrom(chrom.to_owned(), val),
                Some(Err(err)) => StateValue::Error(err),
            },
            StateValue::DiffChrom(c, v) if switch_chrom => StateValue::Value(c, v),
            StateValue::DiffChrom(c, v) => StateValue::DiffChrom(c, v),
            StateValue::Error(e) => StateValue::Error(e),
            StateValue::Done => StateValue::Done,
        };
        // For sanity, if we're switching chromosomes then we should never have an empty value or say we're in a "different" chromosome
        debug_assert!(
            !(switch_chrom
                && matches!(
                    &self.state_value,
                    StateValue::Empty | StateValue::EmptyValue(..) | StateValue::DiffChrom(..)
                )),
        );
    }

    pub(crate) fn load_state_and_take_value(&mut self) -> Option<Result<S::Value, BedValueError>> {
        let state_value = std::mem::replace(&mut self.state_value, StateValue::Empty);
        let ret;
        (self.state_value, ret) = match state_value {
            StateValue::Empty => match self.stream.next() {
                None => (StateValue::Done, None),
                Some(Ok((chrom, val))) => (StateValue::EmptyValue(chrom.to_owned()), Some(Ok(val))),
                Some(Err(err)) => (StateValue::Done, Some(Err(err))),
            },
            StateValue::Value(c, v) => (StateValue::EmptyValue(c), Some(Ok(v))),
            StateValue::EmptyValue(prev_chrom) => match self.stream.next() {
                None => (StateValue::Done, None),
                Some(Ok((chrom, val))) if prev_chrom == chrom => {
                    (StateValue::EmptyValue(prev_chrom), Some(Ok(val)))
                }
                Some(Ok((chrom, val))) => (StateValue::DiffChrom(chrom.to_owned(), val), None),
                Some(Err(err)) => (StateValue::Done, Some(Err(err))),
            },
            StateValue::DiffChrom(c, v) => (StateValue::DiffChrom(c, v), None),
            StateValue::Error(e) => (StateValue::Done, Some(Err(e))),
            StateValue::Done => (StateValue::Done, None),
        };
        // For sanity, we shouldn't have any error or value (for the current chromosome) stored
        debug_assert!(matches!(
            &self.state_value,
            StateValue::Done | StateValue::EmptyValue(..) | StateValue::DiffChrom(..)
        ));
        ret
    }
}

// The separation here between the "current" state and the shared state comes
// from the observation that once we *start* on a chromosome, we can't move on
// to the next until we've exhausted the current. In this *particular*
// implementation, we don't allow parallel iteration of chromsomes. So, the
// state is either needed *here* or in the main struct.
pub struct BedChromData<S: StreamingBedValues> {
    state: Arc<AtomicCell<Option<BedParserState<S>>>>,
    curr_state: Option<BedParserState<S>>,
    pub(crate) done: bool,
}

impl<S: StreamingBedValues> BedChromData<S> {
    pub(crate) fn load_state(&mut self) -> Option<&mut BedParserState<S>> {
        if self.done {
            return None;
        }
        if self.curr_state.is_none() {
            let opt_state = self.state.swap(None);
            if opt_state.is_none() {
                panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
            }
            self.curr_state = opt_state;
        }
        Some(self.curr_state.as_mut().unwrap())
    }
}

impl<S: StreamingBedValues> Drop for BedChromData<S> {
    fn drop(&mut self) {
        if let Some(state) = self.curr_state.take() {
            self.state.swap(Some(state));
        }
    }
}

#[cfg(all(test, features = "write"))]
mod tests {
    use super::*;
    use std::fs::File;
    use std::path::PathBuf;

    use crate::utils::chromvalues::ChromValues;

    #[test]
    fn test_bed_works() {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("small.bed");
        let f = File::open(dir).unwrap();
        let mut bgp = BedParser::from_bed_file(f);
        macro_rules! check_value {
            ($c:ident $chrom:literal) => {
                assert_eq!($c, $chrom);
            };
            (peek next $group:expr, $start:literal $end:literal $rest:expr) => {
                check_value!(peek $group, $start $end $rest);
                check_value!(next $group, $start $end $rest);
            };
            (peek $group:expr, $start:literal $end:literal $rest:expr) => {
                assert_eq!(
                    &BedEntry {
                        start: $start,
                        end: $end,
                        rest: $rest.to_string()
                    },
                    $group.peek().unwrap().unwrap()
                );
            };
            (next $group:expr, $start:literal $end:literal $rest:expr) => {
                assert_eq!(
                    BedEntry {
                        start: $start,
                        end: $end,
                        rest: $rest.to_string()
                    },
                    $group.next().unwrap().unwrap()
                );
            };
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            check_value!(chrom "chr17");
            check_value!(peek group, 1 100 "test1\t0");
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            check_value!(chrom "chr17");
            check_value!(peek group, 1 100 "test1\t0");
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            check_value!(chrom "chr17");
            check_value!(peek next group, 1 100 "test1\t0");
            check_value!(peek next group, 101 200 "test2\t0");
            check_value!(peek next group, 201 300 "test3\t0");
            assert!(group.peek().is_none());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            check_value!(chrom "chr18");
            check_value!(peek next group, 1 100 "test4\t0");
            check_value!(peek next group, 101 200 "test5\t0");
            assert!(group.peek().is_none());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            check_value!(chrom "chr19");
            check_value!(peek next group, 1 100 "test6\t0");
            assert!(group.peek().is_none());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        assert!(matches!(bgp.next_chrom(), None));
    }

    #[test]
    fn test_bedgraph_works() {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("small.bedGraph");
        let f = File::open(dir).unwrap();
        let mut bgp = BedParser::from_bedgraph_file(f);
        macro_rules! check_value {
            ($c:ident $chrom:literal) => {
                assert_eq!($c, $chrom);
            };
            (peek next $group:expr, $start:literal $end:literal) => {
                check_value!(peek $group, $start $end);
                check_value!(next $group, $start $end);
            };
            (peek $group:expr, $start:literal $end:literal) => {
                assert_eq!(
                    &Value {
                        start: $start,
                        end: $end,
                        value: 0.5,
                    },
                    $group.peek().unwrap().unwrap()
                );
            };
            (next $group:expr, $start:literal $end:literal) => {
                assert_eq!(
                    Value {
                        start: $start,
                        end: $end,
                        value: 0.5,
                    },
                    $group.next().unwrap().unwrap()
                );
            };
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            check_value!(chrom "chr17");
            check_value!(peek group, 1 100);
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            check_value!(chrom "chr17");
            check_value!(peek group, 1 100);
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            check_value!(chrom "chr17");
            check_value!(peek next group, 1 100);
            check_value!(peek next group, 101 200);
            check_value!(peek next group, 201 300);
            assert!(group.peek().is_none());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            check_value!(chrom "chr18");
            check_value!(peek next group, 1 100);
            check_value!(peek next group, 101 200);
            assert!(group.peek().is_none());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            check_value!(chrom "chr19");
            check_value!(peek next group, 1 100);
            assert!(group.peek().is_none());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        assert!(matches!(bgp.next_chrom(), None));
    }
}
