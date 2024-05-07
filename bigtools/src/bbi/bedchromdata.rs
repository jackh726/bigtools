//! The types here (`BedParserStreamingIterator` and `BedParserParallelStreamingIterator`)
//! process incoming bed-like data and process into bigWig and bigBed files.
//!
//! `BedParserStreamingIterator` processes the data serially, checking for out
//! of order chromosomes. `BedParserParallelStreamingIterator`, on the other
//! hand, is more complicated wrapper and will queue up to 4 extra chromosomes
//! to be processed concurrently.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

use tokio::runtime::Runtime;

use crate::bed::bedparser::{
    parse_bed, parse_bedgraph, BedFileStream, BedInfallibleIteratorStream, BedIteratorStream,
    BedValueError, Parser, StreamingBedValues,
};
use crate::utils::file_view::FileView;
use crate::utils::streaming_linereader::StreamingLineReader;
use crate::{BedEntry, ChromData, ChromProcess, ProcessChromError, Value};

pub struct BedParserStreamingIterator<S: StreamingBedValues> {
    bed_data: S,
    allow_out_of_order_chroms: bool,
}

impl<S: StreamingBedValues> BedParserStreamingIterator<S> {
    pub fn new(bed_data: S, allow_out_of_order_chroms: bool) -> Self {
        BedParserStreamingIterator {
            bed_data,
            allow_out_of_order_chroms,
        }
    }
}

impl<R: Read> BedParserStreamingIterator<BedFileStream<BedEntry, BufReader<R>>> {
    pub fn from_bed_file(file: R, allow_out_of_order_chroms: bool) -> Self {
        BedParserStreamingIterator::new(
            BedFileStream {
                bed: StreamingLineReader::new(BufReader::new(file)),
                parse: parse_bed,
            },
            allow_out_of_order_chroms,
        )
    }
}

impl<R: Read> BedParserStreamingIterator<BedFileStream<Value, BufReader<R>>> {
    pub fn from_bedgraph_file(file: R, allow_out_of_order_chroms: bool) -> Self {
        BedParserStreamingIterator::new(
            BedFileStream {
                bed: StreamingLineReader::new(BufReader::new(file)),
                parse: parse_bedgraph,
            },
            allow_out_of_order_chroms,
        )
    }
}

impl<
        V: Clone,
        E: Into<BedValueError>,
        C: Into<String> + for<'a> PartialEq<&'a str>,
        I: Iterator<Item = Result<(C, V), E>>,
    > BedParserStreamingIterator<BedIteratorStream<V, I>>
{
    pub fn wrap_iter(iter: I, allow_out_of_order_chroms: bool) -> Self {
        BedParserStreamingIterator::new(
            BedIteratorStream { iter, curr: None },
            allow_out_of_order_chroms,
        )
    }
}

impl<V: Clone, C: Into<String> + for<'a> PartialEq<&'a str>, I: Iterator<Item = (C, V)>>
    BedParserStreamingIterator<BedInfallibleIteratorStream<V, I>>
{
    pub fn wrap_infallible_iter(iter: I, allow_out_of_order_chroms: bool) -> Self {
        BedParserStreamingIterator::new(
            BedInfallibleIteratorStream { iter, curr: None },
            allow_out_of_order_chroms,
        )
    }
}

impl<S: StreamingBedValues> ChromData for BedParserStreamingIterator<S> {
    type Value = S::Value;
    type Error = BedValueError;

    fn process_to_bbi<
        P: ChromProcess<Value = Self::Value>,
        StartProcessing: FnMut(String) -> Result<P, ProcessChromError<Self::Error>>,
        Advance: FnMut(P) -> Result<(), ProcessChromError<Self::Error>>,
    >(
        &mut self,
        runtime: &Runtime,
        start_processing: &mut StartProcessing,
        advance: &mut Advance,
    ) -> Result<(), ProcessChromError<Self::Error>> {
        runtime.block_on(async move {
            let mut state: Option<(String, P, Option<Result<(&str, S::Value), BedValueError>>)> = None;
            loop {
                let (curr_value, new_state) = match state {
                    Some((c, p, Some(v))) => (Some(v), Some((c, p))),
                    Some((c, p, None)) => (self.bed_data.next(), Some((c, p))),
                    None => (self.bed_data.next(), None),
                };
                state = match (new_state, curr_value) {
                    // The next value is an error, but we never started
                    (None, Some(Err(e))) => return Err(ProcessChromError::SourceError(e)),
                    // There are no values at all
                    (None, None) => return Ok(()),
                    // There are no more values
                    (Some(state), None) => {
                        advance(state.1)?;
                        return Ok(());
                    }
                    // The next value is an error and we have seen values before
                    (Some(state), Some(Err(e))) => {
                        // We *can* do anything since we've encountered an error.
                        // We'll go ahead and try to finish what we can, before we return.
                        advance(state.1)?;
                        return Err(ProcessChromError::SourceError(e));
                    }
                    // The next value is the first
                    (None, Some(Ok((chrom, val)))) => {
                        let chrom = chrom.to_string();
                        let mut p = start_processing(chrom.clone())?;
                        let next_val = self.bed_data.next();
                        let next_value = match &next_val {
                            Some(Ok(v)) if v.0 == chrom => Some(&v.1),
                            _ => None,
                        };
                        p.do_process(val, next_value).await?;
                        Some((chrom, p, next_val))
                    }
                    // The next value is the same chromosome
                    (Some((prev_chrom, mut p)), Some(Ok((chrom, val)))) if chrom == &prev_chrom => {
                        let next_val = self.bed_data.next();
                        let next_value = match &next_val {
                            Some(Ok(v)) if v.0 == prev_chrom => Some(&v.1),
                            _ => None,
                        };
                        p.do_process(val, next_value).await?;
                        Some((prev_chrom, p, next_val))
                    }
                    // The next value is a different chromosome
                    (Some((prev_chrom, p)), Some(Ok((chrom, val)))) => {
                        // TODO: test this correctly fails
                        if !self.allow_out_of_order_chroms && prev_chrom.as_str() >= chrom {
                            return Err(ProcessChromError::SourceError(BedValueError::InvalidInput("Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.".to_string())));
                        }
                        advance(p)?;

                        let chrom = chrom.to_string();
                        let mut p = start_processing(chrom.clone())?;
                        let next_val = self.bed_data.next();
                        let next_value = match &next_val {
                            Some(Ok(v)) if v.0 == chrom => Some(&v.1),
                            _ => None,
                        };

                        p.do_process(val, next_value).await?;
                        Some((chrom, p, next_val))
                    }
                };
            }
        })
    }
}

pub struct BedParserParallelStreamingIterator<V> {
    allow_out_of_order_chroms: bool,

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

            chrom_indices,
            parse_fn,
            path,
        }
    }
}

impl<V: Send + 'static> ChromData for BedParserParallelStreamingIterator<V> {
    type Value = V;
    type Error = BedValueError;

    fn process_to_bbi<
        P: ChromProcess<Value = Self::Value> + Send + 'static,
        StartProcessing: FnMut(String) -> Result<P, ProcessChromError<Self::Error>>,
        Advance: FnMut(P) -> Result<(), ProcessChromError<Self::Error>>,
    >(
        &mut self,
        runtime: &Runtime,
        start_processing: &mut StartProcessing,
        advance: &mut Advance,
    ) -> Result<(), ProcessChromError<Self::Error>> {
        let mut remaining = true;
        let mut queued_reads: VecDeque<_> = VecDeque::new();
        loop {
            while remaining && queued_reads.len() < (4 + 1) {
                let (curr, next) = match self.chrom_indices.pop() {
                    Some(c) => (c, self.chrom_indices.last()),
                    None => {
                        remaining = false;
                        break;
                    }
                };
                next.map(|n| assert!(curr.1 != n.1));
                // TODO: test this correctly fails
                if !self.allow_out_of_order_chroms && next.map(|n| curr.1 > n.1).unwrap_or(false) {
                    return Err(ProcessChromError::SourceError(BedValueError::InvalidInput(
                        "Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`."
                            .to_string(),
                    )));
                }

                let file = match File::open(&self.path) {
                    Ok(f) => f,
                    Err(err) => return Err(ProcessChromError::SourceError(err.into())),
                };
                let file = FileView::new(file, curr.0, next.map(|n| n.0).unwrap_or(u64::MAX))?;
                let mut stream = BedFileStream {
                    bed: StreamingLineReader::new(BufReader::new(file)),
                    parse: self.parse_fn,
                };

                let mut p = start_processing(curr.1.clone())?;
                let curr_chrom = curr.1.clone();
                let data: tokio::task::JoinHandle<Result<P, ProcessChromError<BedValueError>>> =
                    runtime.spawn(async move {
                        let mut next_val: Option<Result<(&str, V), BedValueError>> = None;

                        loop {
                            let curr_value = match next_val.take() {
                                Some(v) => Some(v),
                                None => stream.next(),
                            };
                            next_val = match curr_value {
                                // The next value is an error
                                Some(Err(e)) => return Err(ProcessChromError::SourceError(e)),
                                None => return Ok(p),
                                Some(Ok((chrom, _))) if chrom != curr_chrom => {
                                    return Err(ProcessChromError::InvalidInput(
                                        "File is not sorted.".to_string(),
                                    ));
                                }
                                Some(Ok((_, val))) => {
                                    let next_val = stream.next();
                                    let next_value = match &next_val {
                                        Some(Ok(v)) if v.0 == curr_chrom => Some(&v.1),
                                        _ => None,
                                    };
                                    p.do_process(val, next_value).await?;
                                    next_val
                                }
                            };
                        }
                    });
                queued_reads.push_back(data);
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

                Ok(())
            }
        }
        let mut start_processing = |_: String| Ok(TestChromProcess::create(()));
        let mut advance = |p: TestChromProcess| {
            counts.push(p.count);
            let _ = p.destroy();
            Ok(())
        };
        chsi.process_to_bbi(&runtime, &mut start_processing, &mut advance)
            .unwrap();
        assert_eq!(counts, vec![200, 200, 200, 200, 200, 2000]);

        Ok(())
    }
}
