use std::borrow::BorrowMut;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::vec::Vec;

use byteordered::{ByteOrdered, Endianness};
use thiserror::Error;

use crate::bbi::{BBIFile, Summary, Value, ZoomRecord};
use crate::bbiread::{
    get_block_data, read_info, BBIFileInfo, BBIFileReadInfoError, BBIRead, BBIReadError, Block,
    ChromAndSize, ZoomIntervalIter,
};
use crate::utils::reopen::{Reopen, ReopenableFile, SeekableRead};
use crate::{ChromIdNotFound, CirTreeSearchError};

struct IntervalIter<I, R, B>
where
    I: Iterator<Item = Block> + Send,
    R: SeekableRead,
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
    R: SeekableRead,
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

#[derive(Debug, Error)]
pub enum BigWigReadAttachError {
    #[error("NotABigWig")]
    NotABigWig,
    #[error("InvalidChroms")]
    InvalidChroms,
    #[error("{}", .0)]
    IoError(io::Error),
}

impl From<io::Error> for BigWigReadAttachError {
    fn from(error: io::Error) -> Self {
        BigWigReadAttachError::IoError(error)
    }
}

impl From<BBIFileReadInfoError> for BigWigReadAttachError {
    fn from(error: BBIFileReadInfoError) -> Self {
        match error {
            BBIFileReadInfoError::UnknownMagic => BigWigReadAttachError::NotABigWig,
            BBIFileReadInfoError::InvalidChroms => BigWigReadAttachError::InvalidChroms,
            BBIFileReadInfoError::IoError(e) => BigWigReadAttachError::IoError(e),
        }
    }
}

pub struct BigWigRead<R> {
    pub info: BBIFileInfo,
    read: R,
}

impl<R: Reopen> Reopen for BigWigRead<R> {
    fn reopen(&self) -> io::Result<Self> {
        Ok(BigWigRead {
            info: self.info.clone(),
            read: self.read.reopen()?,
        })
    }
}

impl<R: SeekableRead> BBIRead for BigWigRead<R> {
    type Read = R;

    fn get_info(&self) -> &BBIFileInfo {
        &self.info
    }

    fn reader(&mut self) -> &mut R {
        &mut self.read
    }

    fn get_chroms(&self) -> Vec<ChromAndSize> {
        self.info
            .chrom_info
            .iter()
            .map(|c| ChromAndSize {
                name: c.name.clone(),
                length: c.length,
            })
            .collect::<Vec<_>>()
    }
}

impl BigWigRead<ReopenableFile> {
    pub fn from_file_and_attach(path: &str) -> Result<Self, BigWigReadAttachError> {
        let reopen = ReopenableFile {
            path: path.to_string(),
            file: File::open(path)?,
        };
        let b = BigWigRead::from(reopen);
        if b.is_err() {
            eprintln!("Error when opening: {}", path);
        }
        b
    }
}

#[derive(Error, Debug)]
pub enum ZoomIntervalError {
    #[error("The passed reduction level was not found")]
    ReductionLevelNotFound,
    #[error("{}", .0)]
    BBIReadError(BBIReadError),
}

impl From<ChromIdNotFound> for ZoomIntervalError {
    fn from(e: ChromIdNotFound) -> Self {
        ZoomIntervalError::BBIReadError(BBIReadError::InvalidChromosome(e.0))
    }
}

impl From<CirTreeSearchError> for ZoomIntervalError {
    fn from(e: CirTreeSearchError) -> Self {
        ZoomIntervalError::BBIReadError(BBIReadError::CirTreeSearchError(e))
    }
}

impl<R> BigWigRead<R>
where
    R: SeekableRead,
{
    pub fn from(mut read: R) -> Result<Self, BigWigReadAttachError> {
        let info = read_info(&mut read)?;
        match info.filetype {
            BBIFile::BigWig => {}
            _ => return Err(BigWigReadAttachError::NotABigWig),
        }

        Ok(BigWigRead { info, read })
    }

    pub fn inner_read(&self) -> &R {
        &self.read
    }

    /// Does *not* check if the passed `R` matches the provided info (including if the `R` is a bigWig at all!)
    pub fn with_info(info: BBIFileInfo, read: R) -> Self {
        BigWigRead {
            info: info,
            read: read,
        }
    }

    pub fn get_summary(&mut self) -> io::Result<Summary> {
        let endianness = self.info.header.endianness;
        let summary_offset = self.info.header.total_summary_offset;
        let data_offset = self.info.header.full_data_offset;
        let reader = self.reader();
        let mut reader = ByteOrdered::runtime(reader, endianness);
        reader.seek(SeekFrom::Start(summary_offset))?;
        let bases_covered = reader.read_u64()?;
        let min_val = reader.read_f64()?;
        let max_val = reader.read_f64()?;
        let sum = reader.read_f64()?;
        let sum_squares = reader.read_f64()?;
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

    pub fn get_interval<'a>(
        &'a mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> Result<impl Iterator<Item = Result<Value, BBIReadError>> + 'a, BBIReadError> {
        let chrom = self.info.chrom_id(chrom_name)?;
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
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

    pub fn get_interval_move(
        mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> Result<impl Iterator<Item = Result<Value, BBIReadError>>, BBIReadError> {
        let chrom = self.info.chrom_id(chrom_name)?;
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
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

    pub fn get_zoom_interval<'a>(
        &'a mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
        reduction_level: u32,
    ) -> Result<impl Iterator<Item = Result<ZoomRecord, BBIReadError>> + 'a, ZoomIntervalError>
    {
        let chrom = self.info.chrom_id(chrom_name)?;
        let zoom_header = match self
            .info
            .zoom_headers
            .iter()
            .find(|h| h.reduction_level == reduction_level)
        {
            Some(h) => h,
            None => return Err(ZoomIntervalError::ReductionLevelNotFound),
        };

        let index_offset = zoom_header.index_offset;
        let blocks = self.search_cir_tree(index_offset, chrom_name, start, end)?;
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
        let blocks = self
            .get_overlapping_blocks(chrom_name, start, end)
            .map_err(|e| BBIReadError::CirTreeSearchError(e))?;
        let mut values = vec![std::f32::NAN; (end - start) as usize];
        use crate::utils::tell::Tell;
        let mut known_offset = self.reader().tell()?;
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

fn get_block_values<R: SeekableRead>(
    bigwig: &mut BigWigRead<R>,
    block: Block,
    known_offset: &mut u64,
    chrom: u32,
    start: u32,
    end: u32,
) -> Result<Option<std::vec::IntoIter<Value>>, BBIReadError> {
    let mut block_data_mut = get_block_data(bigwig, &block, *known_offset)?;

    use bytes::Buf;
    use bytes::BytesMut;

    let mut bytes_header = BytesMut::zeroed(24);
    block_data_mut.read_exact(&mut bytes_header)?;

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
            let mut bytes = vec![0u8; (item_count as usize) * 12];
            block_data_mut.read_exact(&mut bytes)?;
            assert!(bytes.len() == (item_count as usize) * 12);
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
            let mut bytes = BytesMut::zeroed((item_count as usize) * 8);
            block_data_mut.read_exact(&mut bytes)?;
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
            let mut bytes = BytesMut::zeroed((item_count as usize) * 4);
            block_data_mut.read_exact(&mut bytes)?;
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
