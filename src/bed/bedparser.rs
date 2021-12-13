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

use std::collections::HashMap;
use std::fs::File;
use std::hash::BuildHasher;
use std::io::{self, BufRead, BufReader, Read};
use std::sync::Arc;

use crossbeam_utils::atomic::AtomicCell;

use crate::bigwig::WriteGroupsError;
use crate::bigwig::{BedEntry, Value};
use crate::utils::chromvalues::ChromValues;
use crate::utils::idmap::IdMap;
use crate::utils::streaming_linereader::StreamingLineReader;
use crate::{ChromData, ChromDataState};

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

// FIXME(perf): should this be a generic? Speed test
pub type Parser<V> = Box<dyn for<'a> Fn(&'a str) -> Option<io::Result<(&'a str, V)>> + Send>;

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

#[derive(Debug)]
enum ChromOpt {
    None,
    Same,
    Diff(String),
}

#[derive(Debug)]
struct BedParserState<S: StreamingBedValues> {
    stream: S,
    curr_chrom: Option<String>,
    curr_val: Option<S::Value>,
    next_chrom: ChromOpt,
    next_val: Option<S::Value>,
}

impl<S: StreamingBedValues> BedParser<S> {
    pub fn new(stream: S) -> Self {
        let state = BedParserState {
            stream,
            curr_chrom: None,
            next_chrom: ChromOpt::None,
            curr_val: None,
            next_val: None,
        };
        BedParser {
            state: Arc::new(AtomicCell::new(Some(state))),
        }
    }
}

impl BedParser<BedFileStream<BedEntry, BufReader<File>>> {
    pub fn from_bed_file(file: File) -> Self {
        let parse: Parser<BedEntry> = Box::new(|s: &str| {
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
        });
        BedParser::new(BedFileStream {
            bed: StreamingLineReader::new(BufReader::new(file)),
            parse,
        })
    }
}

impl<R: Read> BedParser<BedFileStream<Value, BufReader<R>>> {
    pub fn from_bedgraph_file(file: R) -> Self {
        let parse: Parser<Value> = Box::new(|s: &str| {
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
        });
        BedParser::new(BedFileStream {
            bed: StreamingLineReader::new(BufReader::new(file)),
            parse,
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
    pub fn next_chrom(&mut self) -> Option<io::Result<(String, BedChromData<S>)>> {
        let mut state = self.state.swap(None).expect("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
        if state.next_val.is_none() {
            match state.advance_state(false) {
                Ok(()) => {}
                Err(e) => return Some(Err(e)),
            }
        }

        let next_chrom = match &state.next_chrom {
            ChromOpt::Diff(real_chrom) => Some(real_chrom),
            ChromOpt::Same => state.curr_chrom.as_ref(),
            ChromOpt::None => None,
        };
        let ret = match next_chrom {
            None => None,
            Some(chrom) => {
                let group = BedChromData {
                    state: self.state.clone(),
                    curr_state: None,
                    done: false,
                };
                Some(Ok((chrom.to_owned(), group)))
            }
        };
        self.state.swap(Some(state));
        ret
    }
}

impl<S: StreamingBedValues> BedParserState<S> {
    fn advance_state(&mut self, replace_current: bool) -> io::Result<()> {
        self.curr_val = self.next_val.take();
        match std::mem::replace(&mut self.next_chrom, ChromOpt::None) {
            ChromOpt::Diff(real_chrom) => {
                self.curr_chrom.replace(real_chrom);
            }
            ChromOpt::Same => {}
            ChromOpt::None => {
                self.curr_chrom = None;
            }
        }

        if let Some(next) = self.stream.next() {
            let (chrom, v) = next?;
            self.next_val.replace(v);
            if let Some(curr_chrom) = &self.curr_chrom {
                if curr_chrom != chrom {
                    self.next_chrom = ChromOpt::Diff(chrom.to_owned());
                } else {
                    self.next_chrom = ChromOpt::Same;
                }
            } else {
                self.next_chrom = ChromOpt::Diff(chrom.to_owned());
            }
        }
        if replace_current && self.curr_val.is_none() && self.next_val.is_some() {
            self.advance_state(false)?;
        }
        Ok(())
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
    type V = S::Value;

    fn next(&mut self) -> Option<io::Result<S::Value>> {
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
        match state.advance_state(true) {
            Ok(()) => {
                if let ChromOpt::Diff(_) = state.next_chrom {
                    self.done = true;
                }
                state.curr_val.take().map(Result::Ok)
            }
            Err(e) => Some(Err(e)),
        }
    }

    fn peek(&mut self) -> Option<&S::Value> {
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
        let state = self.curr_state.as_ref().unwrap();
        state.next_val.as_ref()
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

pub struct BedParserStreamingIterator<S: StreamingBedValues, H: BuildHasher> {
    bed_data: BedParser<S>,
    chrom_map: HashMap<String, u32, H>,
    allow_out_of_order_chroms: bool,
    chrom_ids: Option<IdMap>,
    last_chrom: Option<String>,
}

impl<S: StreamingBedValues, H: BuildHasher> BedParserStreamingIterator<S, H> {
    pub fn new(
        bed_data: BedParser<S>,
        chrom_map: HashMap<String, u32, H>,
        allow_out_of_order_chroms: bool,
    ) -> Self {
        BedParserStreamingIterator {
            bed_data,
            chrom_map,
            allow_out_of_order_chroms,
            chrom_ids: Some(IdMap::default()),
            last_chrom: None,
        }
    }
}

impl<S: StreamingBedValues, H: BuildHasher> ChromData for BedParserStreamingIterator<S, H> {
    type Output = BedChromData<S>;

    fn advance(mut self) -> ChromDataState<Self> {
        match self.bed_data.next_chrom() {
            Some(Err(err)) => ChromDataState::Error(err.into()),
            Some(Ok((chrom, group))) => {
                let chrom_ids = self.chrom_ids.as_mut().unwrap();
                let last = self.last_chrom.replace(chrom.clone());
                if let Some(c) = last {
                    // TODO: test this correctly fails
                    if !self.allow_out_of_order_chroms && c >= chrom {
                        return ChromDataState::Error(WriteGroupsError::InvalidInput("Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.".to_string()));
                    }
                }
                let length = match self.chrom_map.get(&chrom) {
                    Some(length) => *length,
                    None => return ChromDataState::Error(WriteGroupsError::InvalidInput(format!("Input bedGraph contains chromosome that isn't in the input chrom sizes: {}", chrom))),
                };
                let chrom_id = chrom_ids.get_id(&chrom);
                let read_data = (chrom, chrom_id, length, group);

                ChromDataState::Read(read_data, self)
            }
            None => {
                let chrom_ids = self.chrom_ids.take().unwrap();
                ChromDataState::Finished(chrom_ids)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            assert_eq!(chrom, "chr17");
            assert_eq!(
                &BedEntry {
                    start: 1,
                    end: 100,
                    rest: "test1\t0".to_string()
                },
                group.peek().unwrap()
            );
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            assert_eq!(chrom, "chr17");
            assert_eq!(
                &BedEntry {
                    start: 1,
                    end: 100,
                    rest: "test1\t0".to_string()
                },
                group.peek().unwrap()
            );
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            assert_eq!(chrom, "chr17");
            assert_eq!(
                &BedEntry {
                    start: 1,
                    end: 100,
                    rest: "test1\t0".to_string()
                },
                group.peek().unwrap()
            );
            assert_eq!(chrom, "chr17");
            assert_eq!(
                BedEntry {
                    start: 1,
                    end: 100,
                    rest: "test1\t0".to_string()
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(
                &BedEntry {
                    start: 101,
                    end: 200,
                    rest: "test2\t0".to_string()
                },
                group.peek().unwrap()
            );
            assert_eq!(
                &BedEntry {
                    start: 101,
                    end: 200,
                    rest: "test2\t0".to_string()
                },
                group.peek().unwrap()
            );

            assert_eq!(
                BedEntry {
                    start: 101,
                    end: 200,
                    rest: "test2\t0".to_string()
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(
                &BedEntry {
                    start: 201,
                    end: 300,
                    rest: "test3\t0".to_string()
                },
                group.peek().unwrap()
            );

            assert_eq!(
                BedEntry {
                    start: 201,
                    end: 300,
                    rest: "test3\t0".to_string()
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(None, group.peek());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            assert_eq!(chrom, "chr18");
            assert_eq!(
                BedEntry {
                    start: 1,
                    end: 100,
                    rest: "test4\t0".to_string()
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(
                &BedEntry {
                    start: 101,
                    end: 200,
                    rest: "test5\t0".to_string()
                },
                group.peek().unwrap()
            );
            assert_eq!(
                &BedEntry {
                    start: 101,
                    end: 200,
                    rest: "test5\t0".to_string()
                },
                group.peek().unwrap()
            );

            assert_eq!(
                BedEntry {
                    start: 101,
                    end: 200,
                    rest: "test5\t0".to_string()
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(None, group.peek());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            assert_eq!(chrom, "chr19");
            assert_eq!(
                BedEntry {
                    start: 1,
                    end: 100,
                    rest: "test6\t0".to_string()
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(None, group.peek());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        assert!(bgp.next_chrom().is_none());
        Ok(())
    }

    #[test]
    fn test_bedgraph_works() -> io::Result<()> {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("small.bedGraph");
        let f = File::open(dir)?;
        let mut bgp = BedParser::from_bedgraph_file(f);
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            assert_eq!(chrom, "chr17");
            assert_eq!(
                Value {
                    start: 1,
                    end: 100,
                    value: 0.5
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(
                &Value {
                    start: 101,
                    end: 200,
                    value: 0.5
                },
                group.peek().unwrap()
            );
            assert_eq!(
                &Value {
                    start: 101,
                    end: 200,
                    value: 0.5
                },
                group.peek().unwrap()
            );

            assert_eq!(
                Value {
                    start: 101,
                    end: 200,
                    value: 0.5
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(
                &Value {
                    start: 201,
                    end: 300,
                    value: 0.5
                },
                group.peek().unwrap()
            );

            assert_eq!(
                Value {
                    start: 201,
                    end: 300,
                    value: 0.5
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(None, group.peek());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            assert_eq!(chrom, "chr18");
            assert_eq!(
                Value {
                    start: 1,
                    end: 100,
                    value: 0.5
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(
                &Value {
                    start: 101,
                    end: 200,
                    value: 0.5
                },
                group.peek().unwrap()
            );
            assert_eq!(
                &Value {
                    start: 101,
                    end: 200,
                    value: 0.5
                },
                group.peek().unwrap()
            );

            assert_eq!(
                Value {
                    start: 101,
                    end: 200,
                    value: 0.5
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(None, group.peek());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        {
            let (chrom, mut group) = bgp.next_chrom().unwrap().unwrap();
            assert_eq!(chrom, "chr19");
            assert_eq!(
                Value {
                    start: 1,
                    end: 100,
                    value: 0.5
                },
                group.next().unwrap().unwrap()
            );
            assert_eq!(None, group.peek());

            assert!(group.next().is_none());
            assert!(group.peek().is_none());
        }
        assert!(bgp.next_chrom().is_none());
        Ok(())
    }
}
