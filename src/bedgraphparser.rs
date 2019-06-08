use std::fs::File;
use std::io::{self, BufReader};
use std::sync::Arc;

use crate::bigwig::Value;
use crate::streaming_linereader::StreamingLineReader;
use crate::chromvalues::{ChromGroups, ChromValues};

use crossbeam::atomic::AtomicCell;

#[derive(Debug)]
pub struct BedGraphParserState {
    bedgraph: StreamingLineReader<BufReader<File>>,
    curr_chrom: Option<String>,
    curr_val: Option<Value>,
    // None -> No next value
    // Some(None) -> The next chrom is the same as curr_chrom
    // Some(chrom) -> Chrom is different than curr_chrom
    next_chrom: Option<Option<String>>,
    next_val: Option<Value>,
}

impl BedGraphParserState {
    fn advance(&mut self) -> io::Result<()> {
        self.curr_val = self.next_val.take();
        match self.next_chrom.take() {
            Some(next_chrom) => {
                match next_chrom {
                    Some(real_chrom) => {
                        self.curr_chrom.replace(real_chrom);
                    },
                    None => {},
                }
            },
            None => {},
        }

        let l = self.bedgraph.read()?;
        match l {
            Some(line) => {
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
                        self.next_chrom.replace(Some(chrom.to_owned()));
                    } else {
                        self.next_chrom.replace(None);
                    }
                } else {
                    self.next_chrom.replace(Some(chrom.to_owned()));
                }
            },
            None => {
                self.curr_chrom = None;
            },
        }
        if self.curr_val.is_none() && self.next_val.is_some() {
            self.advance()?;
        }
        Ok(())
    }
}

pub struct BedGraphParser {
    state: Arc<AtomicCell<Option<BedGraphParserState>>>,
}

impl BedGraphParser {
    pub fn new(file: File) -> BedGraphParser {
        let bf = BufReader::new(file);
        let state = BedGraphParserState {
            bedgraph: StreamingLineReader::new(bf),
            curr_chrom: None,
            next_chrom: None,
            curr_val: None,
            next_val: None,
        };
        BedGraphParser {
            state: Arc::new(AtomicCell::new(Some(state))),
        }
    }

}

impl ChromGroups<ChromGroup> for BedGraphParser {
    fn next(&mut self) -> io::Result<Option<(String, ChromGroup)>> {
        let opt_state = self.state.swap(None);
        if let None = opt_state {
            panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
        }
        let mut state = opt_state.unwrap();
        if let Some(next_chrom) = state.next_chrom.as_ref() {
            if next_chrom.is_none() {
                panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
            }
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

pub struct ChromGroup {
    state: Arc<AtomicCell<Option<BedGraphParserState>>>,
    curr_state: Option<BedGraphParserState>,
}

impl ChromValues for ChromGroup {
    fn next(&mut self) -> io::Result<Option<Value>> {
        if let None = self.curr_state {
            let opt_state = self.state.swap(None);
            if let None = opt_state {
                panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
            }
            self.curr_state = opt_state;
        }
        let state = self.curr_state.as_mut().unwrap();
        if let Some(val) = state.curr_val.take() {
            return Ok(Some(val));
        }
        if let Some(next_chrom) = state.next_chrom.as_ref() {
            if next_chrom.is_some() {
                return Ok(None);
            }
        }
        state.advance()?;
        Ok(state.curr_val.take())
    }

    fn peek(&mut self) -> Option<&Value> {
        if let None = self.curr_state {
            let opt_state = self.state.swap(None);
            if let None = opt_state {
                panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
            }
            self.curr_state = opt_state;
        }
        let state = self.curr_state.as_ref().unwrap();
        if let Some(next_chrom) = state.next_chrom.as_ref() {
            if next_chrom.is_some() {
                return None;
            }
        }
        return state.next_val.as_ref();
    }
}

impl Drop for ChromGroup {
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
