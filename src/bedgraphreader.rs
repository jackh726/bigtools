use std::sync::Arc;

use crate::bigwig::ValueWithChrom;
use crate::bigwig::Value;

use parking_lot::Mutex;

struct BedGraphReaderState<I> where I: Iterator<Item=ValueWithChrom> {
    inner: Box<I>,
    next_value: Option<Option<ValueWithChrom>>,
}

impl<I> BedGraphReaderState<I> where I: Iterator<Item=ValueWithChrom> {
    fn step(&mut self, chrom: &str) -> Option<ValueWithChrom> {
        let inner = &mut self.inner;
        let next_value = &mut self.next_value;
        let next = next_value.get_or_insert_with(|| inner.next());
        match next {
            None => None,
            Some(next_val) => {
                if next_val.chrom != chrom {
                    None
                } else {
                    next_value.take().unwrap()
                }
            }
        }
    }

    fn step_chrom(&mut self) -> Option<String> {
        let inner = &mut self.inner;
        let next_value = &mut self.next_value;
        let next = next_value.get_or_insert_with(|| inner.next());
        match next {
            None => None,
            Some(next_val) => {
                let chrom = next_val.chrom.clone();
                Some(chrom)
            }
        }
    }
}

pub struct BedGraphReader<I> where I: Iterator<Item=ValueWithChrom> {
    state: Arc<Mutex<BedGraphReaderState<I>>>,
    last_chrom: Option<String>,
}

pub struct BedGraphChromValues<I> where I: Iterator<Item=ValueWithChrom> {
    parent: Arc<Mutex<BedGraphReaderState<I>>>,
    chrom: String,
}

impl<I> BedGraphReader<I> where I: Iterator<Item=ValueWithChrom> {
    pub fn new(iter: I) -> BedGraphReader<I> {
        let state = BedGraphReaderState {
            inner: Box::new(iter),
            next_value: None,
        };
        BedGraphReader {
            state: Arc::new(Mutex::new(state)),
            last_chrom: None,
        }
    }
}

impl<I> Iterator for BedGraphReader<I> where I: Iterator<Item=ValueWithChrom> {
    type Item = (String, BedGraphChromValues<I>);

    fn next(&mut self) -> Option<Self::Item> {
        let mut state = self.state.lock();
        let last_chrom = &mut self.last_chrom;
        let next_chrom = state.step_chrom();
        match next_chrom {
            None => None,
            Some(chrom) => {
                if let Some(last_chrom) = last_chrom {
                    if &chrom == last_chrom {
                        panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
                    }
                }
                let new_chrom = chrom.clone();
                last_chrom.replace(chrom);
                Some((new_chrom.clone(), BedGraphChromValues { chrom: new_chrom, parent: self.state.clone() }))
            }
        }
    }
}

impl<I> Iterator for BedGraphChromValues<I> where I: Iterator<Item=ValueWithChrom> {
    type Item = Value;

    fn next(&mut self) -> Option<Value> {
        let mut parent = self.parent.lock();
        parent.step(&self.chrom).map(|o| Value { start: o.start, end: o.end, value: o.value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate test;

    #[test]
    fn test_works() {
        let vals = vec![
            ValueWithChrom { chrom: String::from("chr1"), start: 0, end: 1, value: 0.5 },
            ValueWithChrom { chrom: String::from("chr1"), start: 1, end: 2, value: 0.3 },
            ValueWithChrom { chrom: String::from("chr2"), start: 0, end: 1, value: 0.5 },
            ValueWithChrom { chrom: String::from("chr2"), start: 1, end: 2, value: 0.3 },
            ValueWithChrom { chrom: String::from("chr3"), start: 0, end: 1, value: 0.5 },
            ValueWithChrom { chrom: String::from("chr3"), start: 1, end: 2, value: 0.3 }
        ];

        let vals_iter = vals.into_iter();

        let reader = BedGraphReader::new(vals_iter);

        for (chrom, values) in reader {
            println!("Next chrom {:?}", chrom);
            for val in values {
                println!("{:?} {:?}", chrom, val);
            }
        }
    }
}
