use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::sync::Arc;

use crate::bigwig::BBIWriteOptions;
use crate::bigwig::BigWigWrite;
use crate::bigwig::ChromGroupRead;
use crate::bigwig::ChromGroupReadStreamingIterator;
use crate::bigwig::Value;
use crate::bigwig::WriteGroupsError;
use crate::idmap::IdMap;
use crate::streaming_linereader::StreamingLineReader;
use crate::chromvalues::{ChromGroups, ChromValues};

use crossbeam::atomic::AtomicCell;

pub fn get_chromgroupstreamingiterator<V: 'static, S: StreamingChromValues + std::marker::Send + 'static>(vals: V, options: BBIWriteOptions, chrom_map: HashMap<String, u32>)
    -> impl ChromGroupReadStreamingIterator
    where V : ChromGroups<ChromGroup<S>> + std::marker::Send {
    struct ChromGroupReadStreamingIteratorImpl<S: StreamingChromValues + std::marker::Send, C: ChromGroups<ChromGroup<S>> + std::marker::Send> {
        chrom_groups: C,
        last_chrom: Option<String>,
        chrom_ids: IdMap<String>,
        pool: futures::executor::ThreadPool,
        options: BBIWriteOptions,
        chrom_map: HashMap<String, u32>,
        _s: std::marker::PhantomData<S>,
    }

    impl<S: StreamingChromValues + std::marker::Send + 'static, C: ChromGroups<ChromGroup<S>> + std::marker::Send> ChromGroupReadStreamingIterator for ChromGroupReadStreamingIteratorImpl<S, C> {
        fn next(&mut self) -> Result<Option<ChromGroupRead>, WriteGroupsError> {
            match self.chrom_groups.next()? {
                Some((chrom, group)) => {
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
                    let chrom_id = self.chrom_ids.get_id(chrom.clone());
                    Ok(Some(BigWigWrite::read_group(chrom, chrom_id, length, group, self.pool.clone(), self.options.clone())?))
                },
                None => Ok(None),
            }
        }
    }

    let group_iter = ChromGroupReadStreamingIteratorImpl {
        chrom_groups: vals,
        last_chrom: None,
        chrom_ids: IdMap::new(),
        pool: futures::executor::ThreadPoolBuilder::new().pool_size(6).create().expect("Unable to create thread pool."),
        options: options.clone(),
        chrom_map: chrom_map,
        _s: std::marker::PhantomData,
    };
    group_iter
}

pub trait StreamingChromValues {
    fn next<'a>(&'a mut self) -> io::Result<Option<(&'a str, u32, u32, f32)>>;
}

pub struct BedGraphStream<B: BufRead> {
    bedgraph: StreamingLineReader<B>
}

impl<B: BufRead> StreamingChromValues for BedGraphStream<B> {
    fn next<'a>(&'a mut self) -> io::Result<Option<(&'a str, u32, u32, f32)>> {
        let l = self.bedgraph.read()?;
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
        let start = split.next().expect("Missing start").parse::<u32>().unwrap();
        let end = split.next().expect("Missing end").parse::<u32>().unwrap();
        let value = split.next().expect("Missing value").parse::<f32>().unwrap();
        Ok(Some((chrom, start, end, value)))
    }
}

pub struct BedGraphIteratorStream<I: Iterator<Item=io::Result<(String, u32, u32, f32)>>> {
    iter: I,
    curr: Option<(String, u32, u32, f32)>,
}

impl<I: Iterator<Item=io::Result<(String, u32, u32, f32)>>> StreamingChromValues for BedGraphIteratorStream<I> {
    fn next<'a>(&'a mut self) -> io::Result<Option<(&'a str, u32, u32, f32)>> {
        use std::ops::Deref;
        self.curr = match self.iter.next() {
            None => return Ok(None),
            Some(v) => Some(v?),
        };
        Ok(self.curr.as_ref().map(|v| (v.0.deref(), v.1, v.2, v.3)))
    }
}

pub struct BedGraphParser<S: StreamingChromValues>{
    state: Arc<AtomicCell<Option<BedGraphParserState<S>>>>,
}

impl<S: StreamingChromValues> BedGraphParser<S> {
    pub fn new(stream: S) -> BedGraphParser<S> {
        let state = BedGraphParserState {
            stream,
            curr_chrom: None,
            next_chrom: ChromOpt::None,
            curr_val: None,
            next_val: None,
        };
        BedGraphParser {
            state: Arc::new(AtomicCell::new(Some(state))),
        }
    }
}

impl BedGraphParser<BedGraphStream<BufReader<File>>> {
    pub fn from_file(file: File) -> BedGraphParser<BedGraphStream<BufReader<File>>> {
        BedGraphParser::new(BedGraphStream { bedgraph: StreamingLineReader::new(BufReader::new(file)) })
    }
}

impl<I: Iterator<Item=io::Result<(String, u32, u32, f32)>>> BedGraphParser<BedGraphIteratorStream<I>> {
    pub fn from_iter(iter: I) -> BedGraphParser<BedGraphIteratorStream<I>> {
        BedGraphParser::new(BedGraphIteratorStream { iter, curr: None })
    }
}

#[derive(Debug)]
enum ChromOpt {
    None,
    Same,
    Diff(String),
}

#[derive(Debug)]
pub struct BedGraphParserState<S: StreamingChromValues> {
    stream: S,
    curr_chrom: Option<String>,
    curr_val: Option<Value>,
    next_chrom: ChromOpt,
    next_val: Option<Value>,
}

impl<S: StreamingChromValues> BedGraphParserState<S> {
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

        if let Some((chrom, start, end, value)) = self.stream.next()? {
            self.next_val.replace(Value { start, end, value });
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

impl<S: StreamingChromValues> ChromGroups<ChromGroup<S>> for BedGraphParser<S> {
    fn next(&mut self) -> io::Result<Option<(String, ChromGroup<S>)>> {
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

pub struct ChromGroup<S: StreamingChromValues> {
    state: Arc<AtomicCell<Option<BedGraphParserState<S>>>>,
    curr_state: Option<BedGraphParserState<S>>,
}

impl<S: StreamingChromValues> ChromValues for ChromGroup<S> {
    fn next(&mut self) -> io::Result<Option<Value>> {
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

    fn peek(&mut self) -> Option<&Value> {
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
        return state.next_val.as_ref();
    }
}

impl<S: StreamingChromValues> Drop for ChromGroup<S> {
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
    extern crate test;

    #[test]
    fn test_works() -> io::Result<()> {
        let f = File::open("/home/hueyj/temp/test.bedGraph")?;
        //let f = File::open("/home/hueyj/temp/final.min.chr17.bedGraph")?;
        let mut bgp = BedGraphParser::from_file(f);
        while let Some((chrom, mut group)) = bgp.next()? {
            println!("Next chrom: {:?}", chrom);
            while let Some(value) = group.next()? {
                println!("Next value: {:?}", value);
                println!("Peek value: {:?}", group.peek());
                println!("Peek value: {:?}", group.peek());
            }
        }
        Ok(())
    }

}
