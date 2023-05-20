use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::vec::Vec;

use byteordered::{ByteOrdered, Endianness};
use thiserror::Error;

use crate::bbiread::{
    get_block_data, read_info, BBIFileInfo, BBIFileReadInfoError, BBIRead, BBIReadError, Block,
    ChromAndSize, MemCachedReader, ZoomIntervalIter,
};
use crate::bigwig::{BBIFile, Summary, Value, ZoomRecord};
use crate::utils::mem_cached_file::{MemCachedRead, CACHE_SIZE};
use crate::utils::seekableread::{Reopen, ReopenableFile, SeekableRead};

struct IntervalIter<'a, I, R, S>
where
    I: Iterator<Item = Block> + Send,
    R: Reopen<S>,
    S: SeekableRead,
{
    bigwig: &'a mut BigWigRead<R, S>,
    known_offset: u64,
    blocks: I,
    // TODO: use type_alias_impl_trait to remove Box
    vals: Option<Box<dyn Iterator<Item = Value> + Send + 'a>>,
    chrom: u32,
    start: u32,
    end: u32,
}

impl<'a, I, R, S> Iterator for IntervalIter<'a, I, R, S>
where
    I: Iterator<Item = Block> + Send,
    R: Reopen<S>,
    S: SeekableRead,
{
    type Item = io::Result<Value>;

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
                        self.bigwig,
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

/// Same as IntervalIter but owned
struct OwnedIntervalIter<I, R, S>
where
    I: Iterator<Item = Block> + Send,
    R: Reopen<S>,
    S: SeekableRead,
{
    bigwig: BigWigRead<R, S>,
    known_offset: u64,
    blocks: I,
    // TODO: use type_alias_impl_trait to remove Box
    vals: Option<Box<dyn Iterator<Item = Value> + Send>>,
    chrom: u32,
    start: u32,
    end: u32,
}

impl<I, R, S> Iterator for OwnedIntervalIter<I, R, S>
where
    I: Iterator<Item = Block> + Send,
    R: Reopen<S>,
    S: SeekableRead,
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
                        &mut self.bigwig,
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
                            return Some(Err(e.into()));
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

pub struct BigWigRead<R, S>
where
    R: Reopen<S>,
    S: SeekableRead,
{
    pub info: BBIFileInfo,
    reopen: R,
    reader: Option<S>,
    cache: HashMap<usize, [u8; CACHE_SIZE]>,
}

impl<R, S> Clone for BigWigRead<R, S>
where
    R: Reopen<S>,
    S: SeekableRead,
{
    fn clone(&self) -> Self {
        BigWigRead {
            info: self.info.clone(),
            reopen: self.reopen.clone(),
            reader: None,
            cache: HashMap::new(),
        }
    }
}

impl<R: Reopen<S>, S: SeekableRead> BBIRead<S> for BigWigRead<R, S> {
    fn get_info(&self) -> &BBIFileInfo {
        &self.info
    }

    fn autosql(&mut self) -> io::Result<String> {
        Ok("".to_string())
    }

    fn ensure_reader(&mut self) -> io::Result<&mut S> {
        if self.reader.is_none() {
            let fp = self.reopen.reopen()?;
            self.reader.replace(fp);
        }
        Ok(self.reader.as_mut().unwrap())
    }

    fn ensure_mem_cached_reader(&mut self) -> io::Result<MemCachedReader<'_, S>> {
        self.ensure_reader()?;
        let inner = self.reader.as_mut().unwrap();
        Ok(MemCachedRead::new(inner, &mut self.cache))
    }

    fn close(&mut self) {
        self.reader.take();
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

impl BigWigRead<ReopenableFile, File> {
    pub fn from_file_and_attach(path: &str) -> Result<Self, BigWigReadAttachError> {
        let reopen = ReopenableFile {
            path: path.to_string(),
        };
        let b = BigWigRead::from(reopen);
        if b.is_err() {
            eprintln!("Error when opening: {}", path);
        }
        b
    }
}

impl<R, S> BigWigRead<R, S>
where
    R: Reopen<S>,
    S: SeekableRead,
{
    pub fn from(reopen: R) -> Result<Self, BigWigReadAttachError> {
        let file = reopen.reopen()?;
        let info = read_info(file)?;
        match info.filetype {
            BBIFile::BigWig => {}
            _ => return Err(BigWigReadAttachError::NotABigWig),
        }

        Ok(BigWigRead {
            info,
            reopen,
            reader: None,
            cache: HashMap::new(),
        })
    }

    pub fn get_summary(&mut self) -> io::Result<Summary> {
        let endianness = self.info.header.endianness;
        let summary_offset = self.info.header.total_summary_offset;
        let data_offset = self.info.header.full_data_offset;
        let reader = self.ensure_reader()?;
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
    ) -> Result<impl Iterator<Item = io::Result<Value>> + Send + 'a, BBIReadError> {
        let chrom = self.info.chrom_id(chrom_name)?;
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
        Ok(IntervalIter {
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
    ) -> Result<impl Iterator<Item = Result<Value, BBIReadError>> + Send, BBIReadError> {
        let chrom = self.info.chrom_id(chrom_name)?;
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
        Ok(OwnedIntervalIter {
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
    ) -> Result<impl Iterator<Item = Result<ZoomRecord, BBIReadError>> + Send + 'a, BBIReadError>
    {
        let chrom = self.info.chrom_id(chrom_name)?;
        let zoom_header = match self
            .info
            .zoom_headers
            .iter()
            .find(|h| h.reduction_level == reduction_level)
        {
            Some(h) => h,
            None => {
                return Err(
                    io::Error::new(io::ErrorKind::Other, "No reduction level found.").into(),
                );
            }
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
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        let mut values = vec![std::f32::NAN; (end - start) as usize];
        use crate::utils::tell::Tell;
        let mut known_offset = self.ensure_reader()?.tell()?;
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

fn get_block_values<R: Reopen<S>, S: SeekableRead>(
    bigwig: &mut BigWigRead<R, S>,
    block: Block,
    known_offset: &mut u64,
    chrom: u32,
    start: u32,
    end: u32,
) -> io::Result<Option<Box<dyn Iterator<Item = Value> + Send>>> {
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
            let mut bytes = BytesMut::zeroed((item_count as usize) * 12);
            block_data_mut.read_exact(&mut bytes)?;
            for _ in 0..item_count {
                // bedgraph
                let (chrom_start, chrom_end, value) = match bigwig.info.header.endianness {
                    Endianness::Big => {
                        let chrom_start = bytes.get_u32();
                        let chrom_end = bytes.get_u32();
                        let value = bytes.get_f32();
                        (chrom_start, chrom_end, value)
                    }
                    Endianness::Little => {
                        let chrom_start = bytes.get_u32_le();
                        let chrom_end = bytes.get_u32_le();
                        let value = bytes.get_f32_le();
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
            return Err(io::Error::new(
                io::ErrorKind::Other,
                format!("Unknown bigwig section type: {}", section_type),
            ));
        }
    }

    *known_offset = block.offset + block.size;
    Ok(Some(Box::new(values.into_iter())))
}
