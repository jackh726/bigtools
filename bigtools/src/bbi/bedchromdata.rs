//! The types here (`BedParserStreamingIterator` and `BedParserParallelStreamingIterator`)
//! ultimately wrap around the a `BedParser` to interface with bigWig and bigBed writing.
//!
//! `BedParserStreamingIterator` is a thin wrapper, which only really has extra checking
//! for out of order chromosomes.
//!
//! `BedParserParallelStreamingIterator` is a more complicated wrapper that will queue up
//! to 4 extra chromosomes to be processed concurrently.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::PathBuf;

use tokio::runtime::{Handle, Runtime};

use crate::bed::bedparser::{
    BedChromData, BedFileStream, BedParser, BedValueError, Parser, StateValue, StreamingBedValues,
};
use crate::utils::chromvalues::ChromValues;
use crate::utils::streaming_linereader::StreamingLineReader;
use crate::{ChromData, ChromProcess, ProcessChromError};

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
    type Values = BedChromData<S>;

    fn process_to_bbi<
        P: ChromProcess<Value = <Self::Values as ChromValues>::Value>,
        StartProcessing: FnMut(String) -> Result<P, ProcessChromError<<Self::Values as ChromValues>::Error>>,
        Advance: FnMut(P) -> Result<(), ProcessChromError<<Self::Values as ChromValues>::Error>>,
    >(
        &mut self,
        runtime: &Runtime,
        start_processing: &mut StartProcessing,
        advance: &mut Advance,
    ) -> Result<(), ProcessChromError<<Self::Values as ChromValues>::Error>> {
        loop {
            match self.bed_data.next_chrom() {
                Some(Ok((chrom, mut group))) => {
                    // First, if we don't want to allow out of order chroms, error here
                    let last = self.last_chrom.replace(chrom.clone());
                    if let Some(c) = last {
                        // TODO: test this correctly fails
                        if !self.allow_out_of_order_chroms && c >= chrom {
                            return Err(ProcessChromError::SourceError(BedValueError::InvalidInput("Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.".to_string())));
                        }
                    }

                    let mut p = start_processing(chrom)?;

                    while let Some(current_val) = group.next() {
                        // If there is a source error, propogate that up
                        let current_val = current_val.map_err(ProcessChromError::SourceError)?;
                        let next_val = match group.peek() {
                            None | Some(Err(_)) => None,
                            Some(Ok(v)) => Some(v),
                        };

                        let read = p.do_process(current_val, next_val);
                        runtime.block_on(read)?;
                    }

                    advance(p)?;
                }
                Some(Err(e)) => return Err(ProcessChromError::SourceError(e)),
                None => break,
            }
        }

        Ok(())
    }
}

pub struct BedParserParallelStreamingIterator<V> {
    allow_out_of_order_chroms: bool,
    last_chrom: Option<String>,

    chrom_indices: Vec<(u64, String)>,
    parse_fn: Parser<V>,
    path: PathBuf,
}

impl<V> BedParserParallelStreamingIterator<V> {
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
        }
    }
}

impl<V: Send + 'static> ChromData for BedParserParallelStreamingIterator<V> {
    type Values = BedChromData<BedFileStream<V, BufReader<File>>>;

    fn process_to_bbi<
        P: ChromProcess<Value = <Self::Values as ChromValues>::Value> + Send + 'static,
        StartProcessing: FnMut(String) -> Result<P, ProcessChromError<<Self::Values as ChromValues>::Error>>,
        Advance: FnMut(P) -> Result<(), ProcessChromError<<Self::Values as ChromValues>::Error>>,
    >(
        &mut self,
        runtime: &Runtime,
        start_processing: &mut StartProcessing,
        advance: &mut Advance,
    ) -> Result<(), ProcessChromError<BedValueError>> {
        let mut remaining = true;
        let mut queued_reads: VecDeque<_> = VecDeque::new();
        loop {
            while remaining && queued_reads.len() < (4 + 1) {
                let curr = match self.chrom_indices.pop() {
                    Some(c) => c,
                    None => {
                        remaining = false;
                        break;
                    }
                };

                let mut file = match File::open(&self.path) {
                    Ok(f) => f,
                    Err(err) => return Err(ProcessChromError::SourceError(err.into())),
                };
                file.seek(SeekFrom::Start(curr.0))?;
                let mut parser = BedParser::new(BedFileStream {
                    bed: StreamingLineReader::new(BufReader::new(file)),
                    parse: self.parse_fn,
                });

                match parser.next_chrom() {
                    Some(Ok((chrom, mut group))) => {
                        let last = self.last_chrom.replace(chrom.clone());
                        if let Some(c) = last {
                            // TODO: test this correctly fails
                            if !self.allow_out_of_order_chroms && c >= chrom {
                                return Err(ProcessChromError::SourceError(BedValueError::InvalidInput("Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.".to_string())));
                            }
                        }

                        let mut p = start_processing(chrom)?;
                        let runtime_handle = runtime.handle().clone();
                        let data: tokio::task::JoinHandle<
                            Result<P, ProcessChromError<BedValueError>>,
                        > = runtime.spawn(async move {
                            while let Some(current_val) = group.next() {
                                // If there is a source error, propogate that up
                                let current_val =
                                    current_val.map_err(ProcessChromError::SourceError)?;
                                let next_val = match group.peek() {
                                    None | Some(Err(_)) => None,
                                    Some(Ok(v)) => Some(v),
                                };

                                let read = p.do_process(current_val, next_val);
                                runtime_handle.block_on(read)?;
                            }
                            Ok(p)
                        });
                        queued_reads.push_back(data);
                    }
                    Some(Err(e)) => return Err(ProcessChromError::SourceError(e)),
                    None => {
                        panic!("Unexpected end of file")
                    }
                }
            }
            let Some(next_chrom) = queued_reads.pop_front() else {
                break;
            };
            let p = runtime.block_on(next_chrom).unwrap()?;
            advance(p)?;
        }

        Ok(())
    }
}

impl<S: StreamingBedValues> ChromValues for BedChromData<S> {
    type Value = S::Value;
    type Error = BedValueError;

    fn next(&mut self) -> Option<Result<Self::Value, Self::Error>> {
        let state = self.load_state()?;
        let ret = state.load_state_and_take_value();
        if matches!(state.state_value, StateValue::DiffChrom(..)) {
            self.done = true;
        }
        ret
    }

    fn peek(&mut self) -> Option<Result<&S::Value, &Self::Error>> {
        let state = self.load_state()?;
        state.load_state(false);
        let ret = match &state.state_value {
            StateValue::Empty => None,
            StateValue::Value(_, val) => Some(Ok(val)),
            StateValue::EmptyValue(_) => None,   // Shouldn't occur
            StateValue::DiffChrom(_, _) => None, // Only `Value` is peekable
            StateValue::Error(err) => Some(Err(err)),
            StateValue::Done => None,
        };
        ret
    }
}

#[cfg(all(test, feature = "write"))]
mod tests {
    use super::*;
    use crate::bed::bedparser::parse_bedgraph;
    use crate::process_internal::ChromProcessCreate;
    use crate::{ProcessChromError, Value};
    use std::fs::File;
    use std::io;
    use std::path::PathBuf;

    #[test]
    fn test_bed_streamingiterator_works() -> io::Result<()> {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("resources/test");
        dir.push("multi_chrom.bedGraph");

        let chrom_indices: Vec<(u64, String)> =
            crate::bed::indexer::index_chroms(File::open(dir.clone())?)?.unwrap();

        let mut chsi = BedParserParallelStreamingIterator::new(
            chrom_indices,
            true,
            PathBuf::from(dir.clone()),
            parse_bedgraph,
        );
        let runtime = tokio::runtime::Builder::new_multi_thread().build().unwrap();
        let mut counts = vec![];
        struct TestChromProcess {
            count: usize,
        }
        impl ChromProcessCreate for TestChromProcess {
            type I = ();
            type Out = ();
            fn create(_: Self::I) -> Self {
                TestChromProcess { count: 0 }
            }
            fn destroy(self) -> Self::Out {}
        }
        impl ChromProcess for TestChromProcess {
            type Value = Value;
            async fn do_process<E: std::error::Error + Send + 'static>(
                &mut self,
                _current_val: Self::Value,
                _next_val: Option<&Self::Value>,
            ) -> Result<(), ProcessChromError<E>> {
                self.count += 1;
                dbg!(self.count);

                Ok(())
            }
        }
        let mut start_processing = |chrom: String| {
            dbg!(chrom);
            Ok(TestChromProcess::create(()))
        };
        let mut advance = |p: TestChromProcess| {
            dbg!(p.count);
            counts.push(p.count);
            let _ = p.destroy();
            Ok(())
        };
        chsi.process_to_bbi(&runtime, &mut start_processing, &mut advance)
            .unwrap();
        assert_eq!(counts, vec![]);

        Ok(())
    }
}
