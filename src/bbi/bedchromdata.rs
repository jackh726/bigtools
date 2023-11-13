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

use crate::bed::bedparser::{
    BedChromData, BedFileStream, BedParser, BedValueError, Parser, StateValue, StreamingBedValues,
};
use crate::utils::chromvalues::ChromValues;
use crate::utils::streaming_linereader::StreamingLineReader;
use crate::{ChromData, ChromDataState, ChromProcessingKey, ProcessChromError};

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

    /// Advancing after `ChromDataState::Finished` has been called will result in a panic.
    fn advance<
        State,
        F: FnMut(
            String,
            BedChromData<S>,
            &mut State,
        ) -> Result<ChromProcessingKey, ProcessChromError<BedValueError>>,
    >(
        &mut self,
        do_read: &mut F,
        state: &mut State,
    ) -> Result<ChromDataState<ChromProcessingKey, BedValueError>, ProcessChromError<BedValueError>>
    {
        Ok(match self.bed_data.next_chrom() {
            Some(Ok((chrom, group))) => {
                // First, if we don't want to allow out of order chroms, error here
                let last = self.last_chrom.replace(chrom.clone());
                if let Some(c) = last {
                    // TODO: test this correctly fails
                    if !self.allow_out_of_order_chroms && c >= chrom {
                        return Ok(ChromDataState::Error(BedValueError::InvalidInput("Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.".to_string())));
                    }
                }

                let read = do_read(chrom, group, state)?;
                ChromDataState::NewChrom(read)
            }
            Some(Err(e)) => ChromDataState::Error(e),
            None => ChromDataState::Finished,
        })
    }
}

pub struct BedParserParallelStreamingIterator<V, E, ChromError> {
    allow_out_of_order_chroms: bool,
    last_chrom: Option<String>,

    chrom_indices: Vec<(u64, String)>,
    parse_fn: Parser<V>,
    path: PathBuf,

    queued_reads: VecDeque<Result<ChromDataState<ChromProcessingKey, ChromError>, E>>,
}

impl<V, E, ChromError> BedParserParallelStreamingIterator<V, E, ChromError> {
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
    for BedParserParallelStreamingIterator<V, ProcessChromError<BedValueError>, BedValueError>
{
    type Values = BedChromData<BedFileStream<V, BufReader<File>>>;

    fn advance<
        State,
        F: FnMut(
            String,
            BedChromData<BedFileStream<V, BufReader<File>>>,
            &mut State,
        ) -> Result<ChromProcessingKey, ProcessChromError<BedValueError>>,
    >(
        &mut self,
        do_read: &mut F,
        state: &mut State,
    ) -> Result<ChromDataState<ChromProcessingKey, BedValueError>, ProcessChromError<BedValueError>>
    {
        let mut begin_next = |_self: &mut Self| -> Result<_, ProcessChromError<BedValueError>> {
            let curr = match _self.chrom_indices.pop() {
                Some(c) => c,
                None => {
                    return Ok(ChromDataState::<_, BedValueError>::Finished);
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
                            return Ok(ChromDataState::Error(BedValueError::InvalidInput("Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.".to_string())));
                        }
                    }

                    let read = do_read(chrom, group, state)?;

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
    use crate::ProcessChromError;
    use std::collections::BTreeMap;
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

        let mut chrom_ids = crate::utils::idmap::IdMap::default();
        let mut key = 0;
        let mut output: BTreeMap<u32, _> = BTreeMap::new();
        let mut do_read = |chrom: String,
                           _: _,
                           output: &mut BTreeMap<u32, _>|
         -> Result<ChromProcessingKey, ProcessChromError<BedValueError>> {
            // Make a new id for the chromosome
            let chrom_id = chrom_ids.get_id(&chrom);

            let curr_key = key;
            key += 1;

            output.insert(curr_key, chrom_id);

            Ok(ChromProcessingKey(curr_key))
        };
        assert!(matches!(
            chsi.advance(&mut do_read, &mut output),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read, &mut output),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read, &mut output),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read, &mut output),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read, &mut output),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read, &mut output),
            Ok(ChromDataState::NewChrom(..))
        ));
        assert!(matches!(
            chsi.advance(&mut do_read, &mut output),
            Ok(ChromDataState::Finished)
        ));

        Ok(())
    }
}
