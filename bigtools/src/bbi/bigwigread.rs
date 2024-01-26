/*!
Provides the interface for reading bigWig files.

## Example
```rust, no_run
# use std::error::Error;
# use std::path::PathBuf;
# use bigtools::BigWigRead;
# use bigtools::BBIRead;
# fn main() -> Result<(), Box<dyn Error>> {
# let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
# dir.push("resources/test");
# let mut bigwig = dir.clone();
# bigwig.push("valid.bigWig");
# let bigwig = bigwig.to_string_lossy();
// First, we open a bigWig using a file name (as a `&str`).
let mut bwread = BigWigRead::open_file(&bigwig)?;

// Then, we could get the chromosomes and lengths
let chroms = bwread.chroms();
assert_eq!(chroms.len(), 1);
assert_eq!(chroms[0].length, 83257441);

// We can get summary data
let summary = bwread.get_summary()?;
assert_eq!(summary.bases_covered, 137894);
assert_eq!(summary.max_val, 14254.0);

// Or we can read data from an interval
let first_interval = bwread
    .get_interval("chr17", 0, 59899)?
    .next()
    .unwrap()
    .unwrap();
assert_eq!(first_interval.start, 59898);
assert_eq!(first_interval.end, 59899);
assert_eq!(first_interval.value, 0.06792);
# Ok(())
# }
```
*/
use std::borrow::BorrowMut;
use std::fs::File;
use std::io::{self, Seek, SeekFrom};
use std::vec::Vec;

use byteordered::{ByteOrdered, Endianness};
use bytes::{Buf, BytesMut};
use thiserror::Error;

use crate::bbi::{BBIFile, Summary, Value, ZoomRecord};
use crate::bbiread::{
    read_info, BBIFileInfo, BBIFileReadInfoError, BBIRead, BBIReadError, Block, ChromInfo,
    ZoomIntervalIter,
};
use crate::internal::BBIReadInternal;
use crate::utils::reopen::{Reopen, ReopenableFile, SeekableRead};
use crate::{search_cir_tree, BBIFileRead, CachedBBIFileRead, ZoomIntervalError};

struct IntervalIter<I, R, B>
where
    I: Iterator<Item = Block> + Send,
    R: BBIFileRead,
    B: BorrowMut<BigWigRead<R>>,
{
    r: std::marker::PhantomData<R>,
    bigwig: B,
    known_offset: u64,
    blocks: I,
    vals: Option<std::vec::IntoIter<Value>>,
    chrom: u32,
    start: u32,
    end: u32,
}

impl<I, R, B> Iterator for IntervalIter<I, R, B>
where
    I: Iterator<Item = Block> + Send,
    R: BBIFileRead,
    B: BorrowMut<BigWigRead<R>>,
{
    type Item = Result<Value, BBIReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match &mut self.vals {
                Some(vals) => match vals.next() {
                    Some(v) => {
                        return Some(Ok(v));
                    }
                    None => {
                        self.vals = None;
                    }
                },
                None => {
                    // TODO: Could minimize this by chunking block reads
                    let current_block = self.blocks.next()?;
                    match get_block_values(
                        self.bigwig.borrow_mut(),
                        current_block,
                        &mut self.known_offset,
                        self.chrom,
                        self.start,
                        self.end,
                    ) {
                        Ok(Some(vals)) => {
                            self.vals = Some(vals);
                        }
                        Ok(None) => {}
                        Err(e) => {
                            return Some(Err(e));
                        }
                    }
                }
            }
        }
    }
}

/// Possible errors encountered when opening a bigWig file to read
#[derive(Debug, Error)]
pub enum BigWigReadOpenError {
    #[error("NotABigWig")]
    NotABigWig,
    #[error("InvalidChroms")]
    InvalidChroms,
    #[error("{}", .0)]
    IoError(io::Error),
}

impl From<io::Error> for BigWigReadOpenError {
    fn from(error: io::Error) -> Self {
        BigWigReadOpenError::IoError(error)
    }
}

impl From<BBIFileReadInfoError> for BigWigReadOpenError {
    fn from(error: BBIFileReadInfoError) -> Self {
        match error {
            BBIFileReadInfoError::UnknownMagic => BigWigReadOpenError::NotABigWig,
            BBIFileReadInfoError::InvalidChroms => BigWigReadOpenError::InvalidChroms,
            BBIFileReadInfoError::IoError(e) => BigWigReadOpenError::IoError(e),
        }
    }
}

/// The struct used to read a bigWig file
pub struct BigWigRead<R> {
    pub(super) info: BBIFileInfo,
    pub(super) read: R,
}

impl<R: Reopen> Reopen for BigWigRead<R> {
    fn reopen(&self) -> io::Result<Self> {
        Ok(BigWigRead {
            info: self.info.clone(),
            read: self.read.reopen()?,
        })
    }
}

impl<R: BBIFileRead> BBIRead for BigWigRead<R> {
    fn info(&self) -> &BBIFileInfo {
        &self.info
    }

    fn chroms(&self) -> &[ChromInfo] {
        &self.info.chrom_info
    }
}

impl<R: BBIFileRead> BBIReadInternal for BigWigRead<R> {
    type Read = R;

    fn reader(&mut self) -> &mut R {
        &mut self.read
    }

    fn reader_and_info(&mut self) -> (&mut Self::Read, &mut BBIFileInfo) {
        (&mut self.read, &mut self.info)
    }
}

impl<R> BigWigRead<R> {
    /// Get basic info about this bigBed
    pub fn info(&self) -> &BBIFileInfo {
        &self.info
    }

    /// Gets the chromosomes present in this bigBed
    pub fn chroms(&self) -> &[ChromInfo] {
        &self.info.chrom_info
    }
}

impl BigWigRead<ReopenableFile> {
    /// Opens a new `BigWigRead` from a given path as a file.
    pub fn open_file(path: &str) -> Result<Self, BigWigReadOpenError> {
        let reopen = ReopenableFile {
            path: path.to_string(),
            file: File::open(path)?,
        };
        let b = BigWigRead::open(reopen);
        if b.is_err() {
            eprintln!("Error when opening: {}", path);
        }
        b
    }
}

impl<R> BigWigRead<R>
where
    R: SeekableRead,
{
    /// Converts this `BigWigRead`` to where the `BBIFileRead` caches index
    /// access and block data
    pub fn cached(self) -> BigWigRead<CachedBBIFileRead<R>> {
        let read = CachedBBIFileRead::new(self.read);
        BigWigRead {
            read,
            info: self.info,
        }
    }
}

impl<R> BigWigRead<R>
where
    R: BBIFileRead,
{
    /// Opens a new `BigWigRead` with for a given type that implements both `Read` and `Seek`
    pub fn open(mut read: R) -> Result<Self, BigWigReadOpenError> {
        let info = read_info(&mut read)?;
        match info.filetype {
            BBIFile::BigWig => {}
            _ => return Err(BigWigReadOpenError::NotABigWig),
        }

        Ok(BigWigRead { info, read })
    }

    /// Does *not* check if the passed `R` matches the provided info (including if the `R` is a bigWig at all!)
    pub fn with_info(info: BBIFileInfo, read: R) -> Self {
        BigWigRead { info, read }
    }

    /// Gets a reference to the inner `R` type, in order to access any info
    pub fn inner_read(&self) -> &R {
        &self.read
    }

    /// Returns the summary data from bigWig
    ///
    /// Note: For version 1 of bigWigs, there is no total summary. In that
    /// case, 0 is returned for all of the summary except total items. If this
    /// matters to you, you can check the version using
    /// `info().header.version > 1`.
    pub fn get_summary(&mut self) -> io::Result<Summary> {
        let endianness = self.info.header.endianness;
        let summary_offset = self.info.header.total_summary_offset;
        let data_offset = self.info.header.full_data_offset;
        let reader = self.reader().raw_reader();
        let mut reader = ByteOrdered::runtime(reader, endianness);
        let (bases_covered, min_val, max_val, sum, sum_squares) = if summary_offset != 0 {
            reader.seek(SeekFrom::Start(summary_offset))?;
            (
                reader.read_u64()?,
                reader.read_f64()?,
                reader.read_f64()?,
                reader.read_f64()?,
                reader.read_f64()?,
            )
        } else {
            (0, 0.0, 0.0, 0.0, 0.0)
        };
        reader.seek(SeekFrom::Start(data_offset))?;
        let total_items = reader.read_u64()?;
        Ok(Summary {
            total_items,
            bases_covered,
            min_val,
            max_val,
            sum,
            sum_squares,
        })
    }

    /// For a given chromosome, start, and end, returns an `Iterator` of the
    /// intersecting `Value`s. The resulting iterator takes a mutable reference
    /// of this `BigWigRead`.
    pub fn get_interval<'a>(
        &'a mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> Result<impl Iterator<Item = Result<Value, BBIReadError>> + 'a, BBIReadError> {
        let chrom = self.info.chrom_id(chrom_name)?;
        let cir_tree = self.full_data_cir_tree()?;
        let blocks = search_cir_tree(&self.info, &mut self.read, cir_tree, chrom_name, start, end)?;
        Ok(IntervalIter {
            r: std::marker::PhantomData,
            bigwig: self,
            known_offset: 0,
            blocks: blocks.into_iter(),
            vals: None,
            chrom,
            start,
            end,
        })
    }

    /// For a given chromosome, start, and end, returns an `Iterator` of the
    /// intersecting `Value`s. The resulting iterator takes this `BigWigRead`
    /// by value.
    pub fn get_interval_move(
        mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> Result<impl Iterator<Item = Result<Value, BBIReadError>>, BBIReadError> {
        let chrom = self.info.chrom_id(chrom_name)?;
        let cir_tree = self.full_data_cir_tree()?;
        let blocks = search_cir_tree(&self.info, &mut self.read, cir_tree, chrom_name, start, end)?;
        Ok(IntervalIter {
            r: std::marker::PhantomData,
            bigwig: self,
            known_offset: 0,
            blocks: blocks.into_iter(),
            vals: None,
            chrom,
            start,
            end,
        })
    }

    /// For a given chromosome, start, and end, returns an `Iterator` of the
    /// intersecting `ZoomRecord`s.
    pub fn get_zoom_interval<'a>(
        &'a mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
        reduction_level: u32,
    ) -> Result<impl Iterator<Item = Result<ZoomRecord, BBIReadError>> + 'a, ZoomIntervalError>
    {
        let cir_tree = self.zoom_cir_tree(reduction_level)?;

        let chrom = self.info.chrom_id(chrom_name)?;

        let blocks = search_cir_tree(&self.info, &mut self.read, cir_tree, chrom_name, start, end)?;

        Ok(ZoomIntervalIter::new(
            self,
            blocks.into_iter(),
            chrom,
            start,
            end,
        ))
    }

    /// Returns the values between `start` and `end` as a `Vec<f32>`. Any
    /// positions with no data in the bigWig will be `std::f32::NAN`.
    pub fn values(
        &mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<f32>, BBIReadError> {
        let chrom = self.info.chrom_id(chrom_name)?;
        let cir_tree = self.full_data_cir_tree()?;
        let blocks = search_cir_tree(&self.info, &mut self.read, cir_tree, chrom_name, start, end)?;
        let mut values = vec![std::f32::NAN; (end - start) as usize];
        use crate::utils::tell::Tell;
        let mut known_offset = self.reader().raw_reader().tell()?;
        for block in blocks {
            let block_values = get_block_values(self, block, &mut known_offset, chrom, start, end)?;
            let block_values = match block_values {
                Some(v) => v,
                None => continue,
            };
            for block_value in block_values {
                let block_value_start = (block_value.start - start) as usize;
                let block_value_end = (block_value.end - start) as usize;
                for i in &mut values[block_value_start..block_value_end] {
                    *i = block_value.value
                }
            }
        }
        Ok(values)
    }
}

fn get_block_values<R: BBIFileRead>(
    bigwig: &mut BigWigRead<R>,
    block: Block,
    known_offset: &mut u64,
    chrom: u32,
    start: u32,
    end: u32,
) -> Result<Option<std::vec::IntoIter<Value>>, BBIReadError> {
    let data = bigwig.read.get_block_data(&bigwig.info, &block)?;
    let mut bytes = BytesMut::with_capacity(data.len());
    bytes.extend_from_slice(&data);

    let mut bytes_header = bytes.split_to(24);

    let (chrom_id, chrom_start, item_step, item_span, section_type, item_count) =
        match bigwig.info.header.endianness {
            Endianness::Big => {
                let chrom_id = bytes_header.get_u32();
                let chrom_start = bytes_header.get_u32();
                let _chrom_end = bytes_header.get_u32();
                let item_step = bytes_header.get_u32();
                let item_span = bytes_header.get_u32();
                let section_type = bytes_header.get_u8();
                let _reserved = bytes_header.get_u8();
                let item_count = bytes_header.get_u16();
                (
                    chrom_id,
                    chrom_start,
                    item_step,
                    item_span,
                    section_type,
                    item_count,
                )
            }
            Endianness::Little => {
                let chrom_id = bytes_header.get_u32_le();
                let chrom_start = bytes_header.get_u32_le();
                let _chrom_end = bytes_header.get_u32_le();
                let item_step = bytes_header.get_u32_le();
                let item_span = bytes_header.get_u32_le();
                let section_type = bytes_header.get_u8();
                let _reserved = bytes_header.get_u8();
                let item_count = bytes_header.get_u16_le();
                (
                    chrom_id,
                    chrom_start,
                    item_step,
                    item_span,
                    section_type,
                    item_count,
                )
            }
        };

    let mut values: Vec<Value> = Vec::with_capacity(item_count as usize);

    if chrom_id != chrom {
        return Ok(None);
    }

    match section_type {
        1 => {
            assert!(bytes.len() >= (item_count as usize) * 12);
            for i in 0..(item_count as usize) {
                let istart = i * 12;
                let block_item_data: &[u8; 12] = bytes[istart..istart + 12].try_into().unwrap();
                // bedgraph
                let (chrom_start, chrom_end, value) = match bigwig.info.header.endianness {
                    Endianness::Big => {
                        let chrom_start = u32::from_be_bytes([
                            block_item_data[0],
                            block_item_data[1],
                            block_item_data[2],
                            block_item_data[3],
                        ]);
                        let chrom_end = u32::from_be_bytes([
                            block_item_data[4],
                            block_item_data[5],
                            block_item_data[6],
                            block_item_data[7],
                        ]);
                        let value = f32::from_be_bytes([
                            block_item_data[8],
                            block_item_data[9],
                            block_item_data[10],
                            block_item_data[11],
                        ]);
                        (chrom_start, chrom_end, value)
                    }
                    Endianness::Little => {
                        let chrom_start = u32::from_le_bytes([
                            block_item_data[0],
                            block_item_data[1],
                            block_item_data[2],
                            block_item_data[3],
                        ]);
                        let chrom_end = u32::from_le_bytes([
                            block_item_data[4],
                            block_item_data[5],
                            block_item_data[6],
                            block_item_data[7],
                        ]);
                        let value = f32::from_le_bytes([
                            block_item_data[8],
                            block_item_data[9],
                            block_item_data[10],
                            block_item_data[11],
                        ]);
                        (chrom_start, chrom_end, value)
                    }
                };
                let mut value = Value {
                    start: chrom_start,
                    end: chrom_end,
                    value,
                };
                if value.end >= start && value.start <= end {
                    value.start = value.start.max(start);
                    value.end = value.end.min(end);
                    values.push(value)
                }
            }
        }
        2 => {
            for _ in 0..item_count {
                // variable step
                let (chrom_start, value) = match bigwig.info.header.endianness {
                    Endianness::Big => {
                        let chrom_start = bytes.get_u32();
                        let value = bytes.get_f32();
                        (chrom_start, value)
                    }
                    Endianness::Little => {
                        let chrom_start = bytes.get_u32_le();
                        let value = bytes.get_f32_le();
                        (chrom_start, value)
                    }
                };
                let chrom_end = chrom_start + item_span;
                let mut value = Value {
                    start: chrom_start,
                    end: chrom_end,
                    value,
                };
                if value.end >= start && value.start <= end {
                    value.start = value.start.max(start);
                    value.end = value.end.min(end);
                    values.push(value)
                }
            }
        }
        3 => {
            let mut curr_start = chrom_start;
            for _ in 0..item_count {
                // fixed step
                let value = match bigwig.info.header.endianness {
                    Endianness::Big => {
                        let value = bytes.get_f32();
                        value
                    }
                    Endianness::Little => {
                        let value = bytes.get_f32_le();
                        value
                    }
                };
                let chrom_start = curr_start;
                curr_start += item_step;
                let chrom_end = chrom_start + item_span;
                let mut value = Value {
                    start: chrom_start,
                    end: chrom_end,
                    value,
                };
                if value.end >= start && value.start <= end {
                    value.start = value.start.max(start);
                    value.end = value.end.min(end);
                    values.push(value)
                }
            }
        }
        _ => {
            return Err(BBIReadError::InvalidFile(format!(
                "Unknown bigwig section type: {}",
                section_type
            )))
        }
    }

    *known_offset = block.offset + block.size;
    Ok(Some(values.into_iter()))
}
