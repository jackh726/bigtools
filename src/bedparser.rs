use std::collections::HashMap;
use std::fs::File;
use std::hash::BuildHasher;
use std::io::{self, BufRead, BufReader};
use std::sync::Arc;

use crossbeam::atomic::AtomicCell;
use futures::future::Either;

use crate::bigwig::{BedEntry, Value};
use crate::bigwig::ChromGroupRead;
use crate::bigwig::ChromGroupReadStreamingIterator;
use crate::bigwig::WriteGroupsError;
use crate::chromvalues::{ChromGroups, ChromValues};
use crate::idmap::IdMap;
use crate::streaming_linereader::StreamingLineReader;


pub type ChromGroupReadFunction<C> = Box<dyn Fn(String, u32, u32, C) -> io::Result<ChromGroupRead> + Send>;

pub struct BedParserChromGroupStreamingIterator<V, C: ChromValues<V> + Send, G: ChromGroups<V, C>, H: BuildHasher> {
    chrom_groups: G,
    callable: ChromGroupReadFunction<C>,
    last_chrom: Option<String>,
    chrom_ids: Option<IdMap>,
    chrom_map: HashMap<String, u32, H>,
    _v: std::marker::PhantomData<V>,
    _s: std::marker::PhantomData<C>,
}

impl<V, C: ChromValues<V> + Send, G: ChromGroups<V, C>, H: BuildHasher> BedParserChromGroupStreamingIterator<V, C, G, H> {
    pub fn new(vals: G, chrom_map: HashMap<String, u32, H>, callable: ChromGroupReadFunction<C>) -> Self{
        BedParserChromGroupStreamingIterator {
            chrom_groups: vals,
            callable,
            last_chrom: None,
            chrom_ids: Some(IdMap::default()),
            chrom_map,
            _v: std::marker::PhantomData,
            _s: std::marker::PhantomData,
        }
    }
}


impl<V, C: ChromValues<V> + Send, G: ChromGroups<V, C>, H: BuildHasher> ChromGroupReadStreamingIterator for BedParserChromGroupStreamingIterator<V, C, G, H> {
    fn next(&mut self) -> Result<Option<Either<ChromGroupRead, (IdMap)>>, WriteGroupsError> {
        match self.chrom_groups.next()? {
            Some((chrom, group)) => {
                let chrom_ids = self.chrom_ids.as_mut().unwrap();
                let last = self.last_chrom.replace(chrom.clone());
                if let Some(c) = last {
                    // TODO: test this correctly fails
                    if c >= chrom {
                        return Err(WriteGroupsError::InvalidInput("Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.".to_string()));
                    }
                }
                let length = match self.chrom_map.get(&chrom) {
                    Some(length) => *length,
                    None => return Err(WriteGroupsError::InvalidInput(format!("Input bedGraph contains chromosome that isn't in the input chrom sizes: {}", chrom))),
                };
                let chrom_id = chrom_ids.get_id(&chrom);
                let group = (self.callable)(chrom, chrom_id, length, group)?;
                Ok(Some(Either::Left(group)))
            },
            None => {
                match self.chrom_ids.take() {
                    Some(chrom_ids) => Ok(Some(Either::Right(chrom_ids))),
                    None => Ok(None),
                }
            }
        }
    }
}


pub trait StreamingChromValues<V> {
    fn next<'a>(&'a mut self) -> io::Result<Option<(&'a str, V)>>;
}

pub struct BedStream<V, B: BufRead> {
    bed: StreamingLineReader<B>,
    parse: Box<dyn Fn(std::str::SplitWhitespace) -> V + Send>,
}

impl<V, B: BufRead> StreamingChromValues<V> for BedStream<V, B> {
    fn next<'a>(&'a mut self) -> io::Result<Option<(&'a str, V)>> {
        let l = self.bed.read()?;
        let line = match l {
            Some(line) => line,
            None => return Ok(None),
        };
        let mut split = line.split_whitespace();
        let chrom = match split.next() {
            Some(chrom) => chrom,
            None => {
                return Ok(None);
            },
        };
        let v = (self.parse)(split);
        Ok(Some((chrom, v)))
    }
}

pub struct BedIteratorStream<V: Clone, I: Iterator<Item=io::Result<(String, V)>>> {
    iter: I,
    curr: Option<(String, V)>,
}

impl<V: Clone, I: Iterator<Item=io::Result<(String, V)>>> StreamingChromValues<V> for BedIteratorStream<V, I> {
    fn next<'a>(&'a mut self) -> io::Result<Option<(&'a str, V)>> {
        use std::ops::Deref;
        self.curr = match self.iter.next() {
            None => return Ok(None),
            Some(v) => Some(v?),
        };
        Ok(self.curr.as_ref().map(|v| (v.0.deref(), v.1.clone())))
    }
}

pub struct BedParser<V, S: StreamingChromValues<V>>{
    state: Arc<AtomicCell<Option<BedParserState<V, S>>>>,
}

impl<V, S: StreamingChromValues<V>> BedParser<V, S> {
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

impl BedParser<BedEntry, BedStream<BedEntry, BufReader<File>>> {
    pub fn from_bed_file(file: File) -> Self {
        let parse = |mut split: std::str::SplitWhitespace<'_> | {
            let start = split.next().expect("Missing start").parse::<u32>().unwrap();
            let end = split.next().expect("Missing end").parse::<u32>().unwrap();
            let rest_strings: Vec<&str> = split.collect();
            let rest = &rest_strings[..].join("\t");
            BedEntry { start, end, rest: rest.to_string() }
        };
        BedParser::new(BedStream { bed: StreamingLineReader::new(BufReader::new(file)), parse: Box::new(parse) })
    }
}

impl BedParser<Value, BedStream<Value, BufReader<File>>> {
    pub fn from_bedgraph_file(file: File) -> Self {
        let parse = |mut split: std::str::SplitWhitespace<'_> | {
            let start = split.next().expect("Missing start").parse::<u32>().unwrap();
            let end = split.next().expect("Missing end").parse::<u32>().unwrap();
            let value = split.next().expect("Missing value").parse::<f32>().unwrap();
            Value { start, end, value }
        };
        BedParser::new(BedStream { bed: StreamingLineReader::new(BufReader::new(file)), parse: Box::new(parse) })
    }
}

impl<V: Clone, I: Iterator<Item=io::Result<(String, V)>>> BedParser<V, BedIteratorStream<V, I>> {
    pub fn from_iter(iter: I) -> Self {
        BedParser::new(BedIteratorStream { iter, curr: None })
    }
}


#[derive(Debug)]
enum ChromOpt {
    None,
    Same,
    Diff(String),
}

#[derive(Debug)]
pub struct BedParserState<V, S: StreamingChromValues<V>> {
    stream: S,
    curr_chrom: Option<String>,
    curr_val: Option<V>,
    next_chrom: ChromOpt,
    next_val: Option<V>,
}

impl<V, S: StreamingChromValues<V>> BedParserState<V, S> {
    fn advance(&mut self) -> io::Result<()> {
        self.curr_val = self.next_val.take();
        match std::mem::replace(&mut self.next_chrom, ChromOpt::None) {
            ChromOpt::Diff(real_chrom) => {
                self.curr_chrom.replace(real_chrom);
            },
            ChromOpt::Same => {},
            ChromOpt::None => {
                self.curr_chrom = None;
            },
        }

        if let Some((chrom, v)) = self.stream.next()? {
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
        if self.curr_val.is_none() && self.next_val.is_some() {
            self.advance()?;
        }
        Ok(())
    }
}

impl<V, S: StreamingChromValues<V>> ChromGroups<V, ChromGroup<V, S>> for BedParser<V, S> {
    fn next(&mut self) -> io::Result<Option<(String, ChromGroup<V, S>)>> {
        let mut state = self.state.swap(None).expect("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
        if let ChromOpt::Same = state.next_chrom {
            panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
        }
        state.advance()?;

        let next_chrom = state.curr_chrom.as_ref();
        let ret = match next_chrom {
            None => Ok(None),
            Some(chrom) => {
                let group = ChromGroup { state: self.state.clone(), curr_state: None };
                Ok(Some((chrom.to_owned(), group)))
            },
        };
        self.state.swap(Some(state));
        ret
    }
}

pub struct ChromGroup<V, S: StreamingChromValues<V>> {
    state: Arc<AtomicCell<Option<BedParserState<V, S>>>>,
    curr_state: Option<BedParserState<V, S>>,
}

impl<V, S: StreamingChromValues<V>> ChromValues<V> for ChromGroup<V, S> {
    fn next(&mut self) -> io::Result<Option<V>> {
        if self.curr_state.is_none() {
            let opt_state = self.state.swap(None);
            if opt_state.is_none() {
                panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
            }
            self.curr_state = opt_state;
        }
        let state = self.curr_state.as_mut().unwrap();
        if let Some(val) = state.curr_val.take() {
            return Ok(Some(val));
        }
        if let ChromOpt::Diff(_) = state.next_chrom {
            return Ok(None);
        }
        state.advance()?;
        Ok(state.curr_val.take())
    }

    fn peek(&mut self) -> Option<&V> {
        if self.curr_state.is_none() {
            let opt_state = self.state.swap(None);
            if opt_state.is_none() {
                panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
            }
            self.curr_state = opt_state;
        }
        let state = self.curr_state.as_ref().unwrap();
        if let ChromOpt::Diff(_) = state.next_chrom {
            return None;
        }
        state.next_val.as_ref()
    }
}

impl<V, S: StreamingChromValues<V>> Drop for ChromGroup<V, S> {
    fn drop(&mut self) {
        if let Some(state) = self.curr_state.take() {
            self.state.swap(Some(state));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io;
    use std::path::PathBuf;
    extern crate test;

    #[test]
    fn test_bed_works() -> io::Result<()> {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("small.bed");
        let f = File::open(dir)?;
        let mut bgp = BedParser::from_bed_file(f);
        {
            let (chrom, mut group) = bgp.next()?.unwrap();
            assert_eq!(chrom, "chr17");
            assert_eq!(BedEntry { start: 1, end: 100, rest: "test1\t0".to_string() }, group.next()?.unwrap());
            assert_eq!(&BedEntry { start: 101, end: 200, rest: "test2\t0".to_string() }, group.peek().unwrap());
            assert_eq!(&BedEntry { start: 101, end: 200, rest: "test2\t0".to_string() }, group.peek().unwrap());

            assert_eq!(BedEntry { start: 101, end: 200, rest: "test2\t0".to_string() }, group.next()?.unwrap());
            assert_eq!(&BedEntry { start: 201, end: 300, rest: "test3\t0".to_string() }, group.peek().unwrap());

            assert_eq!(BedEntry { start: 201, end: 300, rest: "test3\t0".to_string() }, group.next()?.unwrap());
            assert_eq!(None, group.peek());

            assert_eq!(None, group.next()?);
            assert_eq!(None, group.peek());
        }
        {
            let (chrom, mut group) = bgp.next()?.unwrap();
            assert_eq!(chrom, "chr18");
            assert_eq!(BedEntry { start: 1, end: 100, rest: "test4\t0".to_string() }, group.next()?.unwrap());
            assert_eq!(&BedEntry { start: 101, end: 200, rest: "test5\t0".to_string() }, group.peek().unwrap());
            assert_eq!(&BedEntry { start: 101, end: 200, rest: "test5\t0".to_string() }, group.peek().unwrap());

            assert_eq!(BedEntry { start: 101, end: 200, rest: "test5\t0".to_string() }, group.next()?.unwrap());
            assert_eq!(None, group.peek());

            assert_eq!(None, group.next()?);
            assert_eq!(None, group.peek());
        }
        {
            let (chrom, mut group) = bgp.next()?.unwrap();
            assert_eq!(chrom, "chr19");
            assert_eq!(BedEntry { start: 1, end: 100, rest: "test6\t0".to_string() }, group.next()?.unwrap());
            assert_eq!(None, group.peek());

            assert_eq!(None, group.next()?);
            assert_eq!(None, group.peek());
        }
        assert!(bgp.next()?.is_none());
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
            let (chrom, mut group) = bgp.next()?.unwrap();
            assert_eq!(chrom, "chr17");
            assert_eq!(Value { start: 1, end: 100, value: 0.5 }, group.next()?.unwrap());
            assert_eq!(&Value { start: 101, end: 200, value: 0.5 }, group.peek().unwrap());
            assert_eq!(&Value { start: 101, end: 200, value: 0.5 }, group.peek().unwrap());

            assert_eq!(Value { start: 101, end: 200, value: 0.5 }, group.next()?.unwrap());
            assert_eq!(&Value { start: 201, end: 300, value: 0.5 }, group.peek().unwrap());

            assert_eq!(Value { start: 201, end: 300, value: 0.5 }, group.next()?.unwrap());
            assert_eq!(None, group.peek());

            assert_eq!(None, group.next()?);
            assert_eq!(None, group.peek());
        }
        {
            let (chrom, mut group) = bgp.next()?.unwrap();
            assert_eq!(chrom, "chr18");
            assert_eq!(Value { start: 1, end: 100, value: 0.5 }, group.next()?.unwrap());
            assert_eq!(&Value { start: 101, end: 200, value: 0.5 }, group.peek().unwrap());
            assert_eq!(&Value { start: 101, end: 200, value: 0.5 }, group.peek().unwrap());

            assert_eq!(Value { start: 101, end: 200, value: 0.5 }, group.next()?.unwrap());
            assert_eq!(None, group.peek());

            assert_eq!(None, group.next()?);
            assert_eq!(None, group.peek());
        }
        {
            let (chrom, mut group) = bgp.next()?.unwrap();
            assert_eq!(chrom, "chr19");
            assert_eq!(Value { start: 1, end: 100, value: 0.5 }, group.next()?.unwrap());
            assert_eq!(None, group.peek());

            assert_eq!(None, group.next()?);
            assert_eq!(None, group.peek());
        }
        assert!(bgp.next()?.is_none());
        Ok(())
    }

}
