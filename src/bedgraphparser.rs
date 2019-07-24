use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::sync::Arc;

use crate::bigwig::BigWigWrite;
use crate::bigwig::BigWigWriteOptions;
use crate::bigwig::ChromGroupRead;
use crate::bigwig::ChromGroupReadStreamingIterator;
use crate::bigwig::Value;
use crate::bigwig::WriteGroupsError;
use crate::idmap::IdMap;
use crate::streaming_linereader::StreamingLineReader;
use crate::chromvalues::{ChromGroups, ChromValues};

use crossbeam::atomic::AtomicCell;

pub fn get_chromgroupstreamingiterator<V: 'static>(vals: V, options: BigWigWriteOptions) -> impl ChromGroupReadStreamingIterator where V : ChromGroups<ChromGroup<BufReader<File>>> + std::marker::Send {
    struct ChromGroupReadStreamingIteratorImpl<C: ChromGroups<ChromGroup<BufReader<File>>>> {
        chrom_groups: C,
        last_chrom: Option<String>,
        chrom_ids: IdMap<String>,
        pool: futures::executor::ThreadPool,
        options: BigWigWriteOptions,
    }

    impl<C: ChromGroups<ChromGroup<BufReader<File>>>> ChromGroupReadStreamingIterator for ChromGroupReadStreamingIteratorImpl<C> {
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
                    let chrom_id = self.chrom_ids.get_id(chrom.clone());
                    Ok(Some(BigWigWrite::read_group(chrom, chrom_id, group, self.pool.clone(), self.options.clone())?))
                },
                None => Ok(None),
            }
        }
    }

    let group_iter = ChromGroupReadStreamingIteratorImpl {
        chrom_groups: vals,
        last_chrom: None,
        chrom_ids: IdMap::new(),
        pool: futures::executor::ThreadPoolBuilder::new().pool_size(1).create().expect("Unable to create thread pool."),
        options: options.clone(),
    };
    group_iter
}

#[derive(Debug)]
enum ChromOpt {
    None,
    Same,
    Diff(String),
}

#[derive(Debug)]
pub struct BedGraphParserState<B: BufRead> {
    bedgraph: StreamingLineReader<B>,
    curr_chrom: Option<String>,
    curr_val: Option<Value>,
    // None -> No next value
    // Some(None) -> The next chrom is the same as curr_chrom
    // Some(chrom) -> Chrom is different than curr_chrom
    next_chrom: ChromOpt,
    next_val: Option<Value>,
}

impl<B: BufRead> BedGraphParserState<B> {
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

        let l = self.bedgraph.read()?;
        if let Some(line) = l {
            let mut split = line.split_whitespace();
            let chrom = match split.next() {
                Some(chrom) => chrom,
                None => {
                    return Ok(());
                },
            };
            let start = split.next().expect("Missing start").parse::<u32>().unwrap();
            let end = split.next().expect("Missing end").parse::<u32>().unwrap();
            let value = split.next().expect("Missing value").parse::<f32>().unwrap();
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

pub struct BedGraphParser<B: BufRead> {
    state: Arc<AtomicCell<Option<BedGraphParserState<B>>>>,
}

impl BedGraphParser<BufReader<File>> {
    pub fn new(file: File) -> BedGraphParser<BufReader<File>> {
        let bf = BufReader::new(file);
        let state = BedGraphParserState {
            bedgraph: StreamingLineReader::new(bf),
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

impl<B: BufRead> ChromGroups<ChromGroup<B>> for BedGraphParser<B> {
    fn next(&mut self) -> io::Result<Option<(String, ChromGroup<B>)>> {
        let opt_state = self.state.swap(None);
        if opt_state.is_none() {
            panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
        }
        let mut state = opt_state.unwrap();
        if let ChromOpt::Same = state.next_chrom {
            panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
        }
        state.advance()?;

        let next_chrom = state.curr_chrom.as_ref();
        let ret = match next_chrom {
            None => Ok(None),
            Some(chrom) => Ok(Some((chrom.to_owned(), ChromGroup { state: self.state.clone(), curr_state: None } ))),

        };
        self.state.swap(Some(state));
        ret
    }
}

pub struct ChromGroup<B: BufRead> {
    state: Arc<AtomicCell<Option<BedGraphParserState<B>>>>,
    curr_state: Option<BedGraphParserState<B>>,
}

impl<B: BufRead> ChromValues for ChromGroup<B> {
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

impl<B: BufRead> Drop for ChromGroup<B> {
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
        let mut bgp = BedGraphParser::new(f);
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
