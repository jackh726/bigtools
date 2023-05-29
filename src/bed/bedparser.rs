//! There are roughly three layers of abstraction here, each with a different purpose.
//!
//! The first layer of abstraction is enscapsulated in the `StreamingBedValues` trait. Briefly,
//! implementors of this trait return "raw" bed-like data. This is the chromosome (as a `&str`) and
//! data specific to each type of bed
//!
//! The second layer of abstraction manages the state information for when values switch from one
//! chromosome to another. The is important because bigwig/bigbed writing is "chunked" by chromosome.
//!
//! The final layer of abstraction is a thin wrapper around the previous to provide some optional
//! error checking and to keep track of the chromosomes seen.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::PathBuf;
use std::sync::Arc;

use crossbeam_utils::atomic::AtomicCell;
use thiserror::Error;

use crate::bigwig::{BedEntry, Value};
use crate::utils::chromvalues::ChromValues;
use crate::utils::streaming_linereader::StreamingLineReader;
use crate::{ChromData, ChromDataState, ChromProcessingFnOutput};

// FIXME: replace with LendingIterator when GATs are thing
/// Essentially a combined lending iterator over the chrom (&str) and remaining
/// values of bed-like data
pub trait StreamingBedValues {
    type Value;

    fn next(&mut self) -> Option<io::Result<(&str, Self::Value)>>;
}

// ---------------
// Bed-like stream
// ---------------

pub type Parser<V> = for<'a> fn(&'a str) -> Option<io::Result<(&'a str, V)>>;

/// Parses a bed-like file
pub struct BedFileStream<V, B> {
    bed: StreamingLineReader<B>,
    parse: Parser<V>,
}

impl<V, B: BufRead> StreamingBedValues for BedFileStream<V, B> {
    type Value = V;

    fn next(&mut self) -> Option<io::Result<(&str, Self::Value)>> {
        let line = match self.bed.read()? {
            Ok(line) => line.trim_end(),
            Err(e) => return Some(Err(e)),
        };
        (self.parse)(line)
    }
}

// Wraps a bed-like Iterator
pub struct BedIteratorStream<V, I> {
    iter: I,
    curr: Option<(String, V)>,
}

impl<V: Clone, I: Iterator<Item = io::Result<(String, V)>>> StreamingBedValues
    for BedIteratorStream<V, I>
{
    type Value = V;

    fn next(&mut self) -> Option<io::Result<(&str, V)>> {
        use std::ops::Deref;
        self.curr = match self.iter.next()? {
            Err(e) => return Some(Err(e)),
            Ok(v) => Some(v),
        };
        self.curr.as_ref().map(|v| Ok((v.0.deref(), v.1.clone())))
    }
}

// ----------------
// State-management
// ----------------

/// A wrapper for "bed-like" data
pub struct BedParser<S: StreamingBedValues> {
    state: Arc<AtomicCell<Option<BedParserState<S>>>>,
}

/// Defines the internal states of bed parsing
enum StateValue<V> {
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
    Error(BedParseError),
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
    fn take_error(&mut self) -> Option<BedParseError> {
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

// Example order of state transitions
// 1) active_chrom: None, next_val: None (creation)
// 2) active_chrom: Some(X), next_val: Some((.., Same)) (load value)
// 3) active_chrom: Some(X), next_val: None (value taken)
// (cycle between for 2 and 3 for all values of a chromosome)
// 4) active_chrom: None, next_val: Some((.., Diff(Y))) (switch chromosome)
// 5) active_chrom: Some(Y), next_val: Some((.. Same)) (load value)
// 6) active_chrom: Some(Y), next_val: None (value taken)
// (cycle between 5 and 6 for all values of a chromosome)
#[derive(Debug)]
struct BedParserState<S: StreamingBedValues> {
    stream: S,
    state_value: StateValue<S::Value>,
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

pub fn parse_bed<'a>(s: &'a str) -> Option<io::Result<(&'a str, BedEntry)>> {
    let mut split = s.splitn(4, '\t');
    let chrom = match split.next() {
        Some(chrom) => chrom,
        None => return None,
    };
    let res = (|| {
        let s = split.next().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Missing start: {:}", s))
        })?;
        let start = s.parse::<u32>().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid start: {:}", s))
        })?;
        let s = split.next().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Missing end: {:}", s))
        })?;
        let end = s.parse::<u32>().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid end: {:}", s))
        })?;
        let rest = split.next().unwrap_or("").to_string();
        Ok((start, end, rest))
    })();
    match res {
        Err(e) => Some(Err(e)),
        Ok((start, end, rest)) => Some(Ok((chrom, BedEntry { start, end, rest }))),
    }
}

pub fn parse_bedgraph<'a>(s: &'a str) -> Option<io::Result<(&'a str, Value)>> {
    let mut split = s.splitn(5, '\t');
    let chrom = match split.next() {
        Some(chrom) => chrom,
        None => return None,
    };
    let res = (|| {
        let s = split.next().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Missing start: {:}", s))
        })?;
        let start = s.parse::<u32>().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid start: {:}", s))
        })?;
        let s = split.next().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Missing end: {:}", s))
        })?;
        let end = s.parse::<u32>().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid end: {:}", s))
        })?;
        let s = split.next().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Missing value: {:}", s))
        })?;
        let value = s.parse::<f32>().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid value: {:}", s))
        })?;
        Ok((start, end, value))
    })();
    match res {
        Err(e) => Some(Err(e)),
        Ok((start, end, value)) => Some(Ok((chrom, Value { start, end, value }))),
    }
}

impl BedParser<BedFileStream<BedEntry, BufReader<File>>> {
    pub fn from_bed_file(file: File) -> Self {
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

impl<V: Clone, I: Iterator<Item = io::Result<(String, V)>>> BedParser<BedIteratorStream<V, I>> {
    pub fn wrap_iter(iter: I) -> Self {
        BedParser::new(BedIteratorStream { iter, curr: None })
    }
}

impl<S: StreamingBedValues> BedParser<S> {
    // This is *valid* to call multiple times for the same chromosome (assuming the
    // `BedChromData` has been dropped), since calling this function doesn't
    // actually advance the state (it will only set `next_val` if it currently is none).
    pub fn next_chrom(&mut self) -> Option<Result<(String, BedChromData<S>), BedParseError>> {
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
    fn load_state(&mut self, switch_chrom: bool) {
        let state_value = std::mem::replace(&mut self.state_value, StateValue::Empty);
        self.state_value = match state_value {
            StateValue::Empty => match self.stream.next() {
                None => StateValue::Done,
                Some(Ok((chrom, val))) => StateValue::Value(chrom.to_owned(), val),
                Some(Err(err)) => StateValue::Error(BedParseError::IoError(err)),
            },
            StateValue::Value(c, v) => StateValue::Value(c, v),
            StateValue::EmptyValue(prev_chrom) => match self.stream.next() {
                None => StateValue::Done,
                Some(Ok((chrom, val))) if switch_chrom || prev_chrom == chrom => {
                    StateValue::Value(prev_chrom, val)
                }
                Some(Ok((chrom, val))) => StateValue::DiffChrom(chrom.to_owned(), val),
                Some(Err(err)) => StateValue::Error(BedParseError::IoError(err)),
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

    fn load_state_and_take_value(&mut self) -> Option<Result<S::Value, BedParseError>> {
        let state_value = std::mem::replace(&mut self.state_value, StateValue::Empty);
        let ret;
        (self.state_value, ret) = match state_value {
            StateValue::Empty => match self.stream.next() {
                None => (StateValue::Done, None),
                Some(Ok((chrom, val))) => (StateValue::EmptyValue(chrom.to_owned()), Some(Ok(val))),
                Some(Err(err)) => (StateValue::Done, Some(Err(BedParseError::IoError(err)))),
            },
            StateValue::Value(c, v) => (StateValue::EmptyValue(c), Some(Ok(v))),
            StateValue::EmptyValue(prev_chrom) => match self.stream.next() {
                None => (StateValue::Done, None),
                Some(Ok((chrom, val))) if prev_chrom == chrom => {
                    (StateValue::EmptyValue(prev_chrom), Some(Ok(val)))
                }
                Some(Ok((chrom, val))) => (StateValue::DiffChrom(chrom.to_owned(), val), None),
                Some(Err(err)) => (StateValue::Done, Some(Err(BedParseError::IoError(err)))),
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

#[derive(Debug, Error)]
pub enum BedParseError {
    #[error("{}", .0)]
    InvalidInput(String),
    #[error("{}", .0)]
    IoError(io::Error),
}

impl From<io::Error> for BedParseError {
    fn from(e: io::Error) -> Self {
        Self::IoError(e)
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
    done: bool,
}

impl<S: StreamingBedValues> ChromValues for BedChromData<S> {
    type Value = S::Value;
    type Error = BedParseError;

    fn next(&mut self) -> Option<Result<Self::Value, Self::Error>> {
        if self.curr_state.is_none() {
            let opt_state = self.state.swap(None);
            if opt_state.is_none() {
                panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
            }
            self.curr_state = opt_state;
        }
        if self.done {
            return None;
        }
        let state = self.curr_state.as_mut().unwrap();
        let ret = state.load_state_and_take_value();
        if matches!(state.state_value, StateValue::DiffChrom(..)) {
            self.done = true;
        }
        ret
    }

    fn peek(&mut self) -> Option<Result<&S::Value, &Self::Error>> {
        if self.curr_state.is_none() {
            let opt_state = self.state.swap(None);
            if opt_state.is_none() {
                panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
            }
            self.curr_state = opt_state;
        }
        if self.done {
            return None;
        }
        let state = self.curr_state.as_mut().unwrap();
        state.load_state(false);
        let ret = match &state.state_value {
            StateValue::Empty => None,
            StateValue::Value(_, val) => Some(Ok(val)),
            StateValue::EmptyValue(_) => None,   // Shouldn't occur
            StateValue::DiffChrom(_, _) => None, // Only `Value` is peekable
            StateValue::Error(err) => Some(Err(err)),
            StateValue::Done => None,
        };
        if matches!(&state.state_value, StateValue::DiffChrom(..)) {
            self.done = true;
        }
        ret
    }
}

impl<S: StreamingBedValues> Drop for BedChromData<S> {
    fn drop(&mut self) {
        if let Some(state) = self.curr_state.take() {
            self.state.swap(Some(state));
        }
    }
}

// ------------------------------------------------
// Chromosome tracking and optional error reporting
// ------------------------------------------------

pub struct BedParserStreamingIterator<S: StreamingBedValues> {
    bed_data: BedParser<S>,
    allow_out_of_order_chroms: bool,
    last_chrom: Option<String>,
}

impl<S: StreamingBedValues> BedParserStreamingIterator<S> {
    pub fn new(bed_data: BedParser<S>, allow_out_of_order_chroms: bool) -> Self {
        BedParserStreamingIterator {
            bed_data,
            allow_out_of_order_chroms,
            last_chrom: None,
        }
    }
}

impl<S: StreamingBedValues> ChromData for BedParserStreamingIterator<S> {
    type Output = BedChromData<S>;

    /// Advancing after `ChromDataState::Finished` has been called will result in a panic.
    fn advance<
        F: FnMut(String, Self::Output) -> io::Result<ChromProcessingFnOutput<Self::Output>>,
    >(
        &mut self,
        do_read: &mut F,
    ) -> io::Result<ChromDataState<Self::Output>> {
        Ok(match self.bed_data.next_chrom() {
            Some(Ok((chrom, group))) => {
                // First, if we don't want to allow out of order chroms, error here
                let last = self.last_chrom.replace(chrom.clone());
                if let Some(c) = last {
                    // TODO: test this correctly fails
                    if !self.allow_out_of_order_chroms && c >= chrom {
                        return Ok(ChromDataState::Error(BedParseError::InvalidInput("Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.".to_string())));
                    }
                }

                let read = do_read(chrom, group)?;
                ChromDataState::NewChrom(read)
            }
            Some(Err(e)) => ChromDataState::Error(e),
            None => ChromDataState::Finished,
        })
    }
}

pub struct BedParserParallelStreamingIterator<V, O: ChromValues> {
    allow_out_of_order_chroms: bool,
    last_chrom: Option<String>,

    chrom_indices: Vec<(u64, String)>,
    parse_fn: Parser<V>,
    path: PathBuf,

    queued_reads: VecDeque<io::Result<ChromDataState<O>>>,
}

impl<V, O: ChromValues> BedParserParallelStreamingIterator<V, O> {
    pub fn new(
        mut chrom_indices: Vec<(u64, String)>,
        allow_out_of_order_chroms: bool,
        path: PathBuf,
        parse_fn: Parser<V>,
    ) -> Self {
        // For speed, we `pop` and go in reverse order. We want forward order,
        // so reverse here.
        chrom_indices.reverse();

        BedParserParallelStreamingIterator {
            allow_out_of_order_chroms,
            last_chrom: None,

            chrom_indices,
            parse_fn,
            path,

            queued_reads: VecDeque::new(),
        }
    }
}

impl<V> ChromData
    for BedParserParallelStreamingIterator<V, BedChromData<BedFileStream<V, BufReader<File>>>>
{
    type Output = BedChromData<BedFileStream<V, BufReader<File>>>;

    fn advance<
        F: FnMut(String, Self::Output) -> io::Result<ChromProcessingFnOutput<Self::Output>>,
    >(
        &mut self,
        do_read: &mut F,
    ) -> io::Result<ChromDataState<Self::Output>> {
        let mut begin_next = |_self: &mut Self| -> io::Result<_> {
            let curr = match _self.chrom_indices.pop() {
                Some(c) => c,
                None => {
                    return Ok(ChromDataState::<Self::Output>::Finished);
                }
            };

            let mut file = match File::open(&_self.path) {
                Ok(f) => f,
                Err(err) => return Ok(ChromDataState::Error(err.into())),
            };
            file.seek(SeekFrom::Start(curr.0))?;
            let mut parser = BedParser::new(BedFileStream {
                bed: StreamingLineReader::new(BufReader::new(file)),
                parse: _self.parse_fn,
            });

            Ok(match parser.next_chrom() {
                Some(Ok((chrom, group))) => {
                    let last = _self.last_chrom.replace(chrom.clone());
                    if let Some(c) = last {
                        // TODO: test this correctly fails
                        if !_self.allow_out_of_order_chroms && c >= chrom {
                            return Ok(ChromDataState::Error(BedParseError::InvalidInput("Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.".to_string())));
                        }
                    }

                    let read = do_read(chrom, group)?;

                    ChromDataState::NewChrom(read)
                }
                Some(Err(e)) => ChromDataState::Error(e),
                None => {
                    panic!("Unexpected end of file")
                }
            })
        };

        while self.queued_reads.len() < (4 + 1)
            && matches!(
                self.queued_reads.back(),
                None | Some(Ok(ChromDataState::NewChrom(..)))
            )
        {
            let next = begin_next(self);
            self.queued_reads.push_back(next);
        }
        self.queued_reads.pop_front().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BBIWriteOptions;
    use std::fs::File;
    use std::io;
    use std::path::PathBuf;

    #[test]
    fn test_bed_works() -> io::Result<()> {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("small.bed");
        let f = File::open(dir)?;
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
        Ok(())
    }

    #[test]
    fn test_bedgraph_works() -> io::Result<()> {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("small.bedGraph");
        let f = File::open(dir)?;
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
        Ok(())
    }

    #[test]
    fn test_bed_streamingiterator_works() -> io::Result<()> {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("multi_chrom.bedGraph");

        let chrom_map = std::collections::HashMap::from([
            ("chr1".to_owned(), 100000),
            ("chr2".to_owned(), 100000),
            ("chr3".to_owned(), 100000),
            ("chr4".to_owned(), 100000),
            ("chr5".to_owned(), 100000),
            ("chr6".to_owned(), 100000),
        ]);

        let chrom_indices: Vec<(u64, String)> =
            crate::bed::indexer::index_chroms(File::open(dir.clone())?)?;

        let mut chsi = BedParserParallelStreamingIterator::new(
            chrom_indices,
            true,
            PathBuf::from(dir.clone()),
            parse_bedgraph,
        );

        let pool = futures::executor::ThreadPoolBuilder::new()
            .pool_size(1)
            .create()
            .expect("Unable to create thread pool.");
        let options = BBIWriteOptions::default();

        let mut chrom_ids = crate::utils::idmap::IdMap::default();
        let mut do_read = |chrom: String,
                           data|
         -> io::Result<
            ChromProcessingFnOutput<BedChromData<BedFileStream<Value, BufReader<File>>>>,
        > {
            let length = match chrom_map.get(&chrom) {
                Some(length) => *length,
                None => return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Input bedGraph contains chromosome that isn't in the input chrom sizes: {}", 0),
                )),
            };
            // Make a new id for the chromosome
            let chrom_id = chrom_ids.get_id(&chrom);

            crate::BigWigWrite::begin_processing_chrom(
                chrom,
                data,
                pool.clone(),
                options,
                chrom_id,
                length,
            )
        };
        assert!(matches!(
            chsi.advance(&mut do_read),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read),
            Ok(ChromDataState::Finished)
        ));

        Ok(())
    }
}
