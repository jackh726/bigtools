use std::sync::Arc;

use crate::bigwig::ValueWithChrom;
use crate::bigwig::Value;

use crossbeam::atomic::AtomicCell;

struct BedGraphReaderState {
    inner: Box<Iterator<Item=ValueWithChrom> + std::marker::Send>,
    next_value: Option<Option<ValueWithChrom>>,
}

impl BedGraphReaderState {
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

pub struct BedGraphReader {
    state: Arc<AtomicCell<Option<BedGraphReaderState>>>,
    last_chrom: Option<String>,
}

pub struct BedGraphChromValues {
    state: Arc<AtomicCell<Option<BedGraphReaderState>>>,
    curr_state: Option<BedGraphReaderState>,
    chrom: String,
}

impl BedGraphReader {
    pub fn new<I: Iterator<Item=ValueWithChrom> + std::marker::Send + 'static>(iter: I) -> BedGraphReader {
        let state = BedGraphReaderState {
            inner: Box::new(iter),
            next_value: None,
        };
        BedGraphReader {
            state: Arc::new(AtomicCell::new(Some(state))),
            last_chrom: None,
        }
    }
}

impl Iterator for BedGraphReader {
    type Item = (String, BedGraphChromValues);

    fn next(&mut self) -> Option<Self::Item> {
        let opt_state = self.state.swap(None);
        if let None = opt_state {
            panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
        }
        let mut state = opt_state.unwrap();
        let last_chrom = &mut self.last_chrom;
        let next_chrom = state.step_chrom();
        let ret = match next_chrom {
            None => None,
            Some(chrom) => {
                if let Some(last_chrom) = last_chrom {
                    if &chrom == last_chrom {
                        panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
                    }
                }
                let new_chrom = chrom.clone();
                last_chrom.replace(chrom);
                Some((new_chrom.clone(), BedGraphChromValues { chrom: new_chrom, state: self.state.clone(), curr_state: None, }))
            }
        };
        self.state.swap(Some(state));
        ret
    }
}

impl Iterator for BedGraphChromValues {
    type Item = Value;

    fn next(&mut self) -> Option<Value> {
        if let None = self.curr_state {
            let opt_state = self.state.swap(None);
            if let None = opt_state {
                panic!("Invalid usage. This iterator does not buffer and all values should be exhausted for a chrom before next() is called.");
            }
            self.curr_state = opt_state;
        }
        let state = self.curr_state.as_mut().unwrap();
        let ret = state.step(&self.chrom).map(|o| Value { start: o.start, end: o.end, value: o.value });
        ret
    }
}

impl Drop for BedGraphChromValues {
    fn drop(&mut self) {
        if let Some(state) = self.curr_state.take() {
            self.state.swap(Some(state));
        }
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

    #[bench]
    fn bench_works(b: &mut test::Bencher) {
        b.iter(|| {
            let vals = vec![
                ValueWithChrom { chrom: String::from("chr1"), start: 0, end: 1, value: 0.5 },
                ValueWithChrom { chrom: String::from("chr1"), start: 1, end: 2, value: 0.3 },
                ValueWithChrom { chrom: String::from("chr1"), start: 3, end: 4, value: 0.5 },
                ValueWithChrom { chrom: String::from("chr1"), start: 4, end: 5, value: 0.3 },
                ValueWithChrom { chrom: String::from("chr3"), start: 0, end: 1, value: 0.5 },
                ValueWithChrom { chrom: String::from("chr3"), start: 1, end: 2, value: 0.3 }
            ];

            let vals_iter = vals.into_iter();

            let reader = BedGraphReader::new(vals_iter);

            for (_chrom, values) in reader {
                for val in values {
                    drop(val);
                }
            }
        })
    }
}
