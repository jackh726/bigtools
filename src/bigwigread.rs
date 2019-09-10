use std::io::{self, Seek, SeekFrom};
use std::io::{BufReader};
use std::fs::File;
use std::vec::Vec;

use byteordered::{ByteOrdered, Endianness};

use crate::seekableread::{Reopen, SeekableRead, ReopenableFile};
use crate::bbiread::{BBIRead, BBIFileReadInfoError, BBIFileInfo, Block, ChromAndSize, ZoomIntervalIter, read_info, get_block_data};
use crate::bigwig::{BBIFile, Value, ZoomRecord};


struct IntervalIter<'a, I, R, S> where I: Iterator<Item=Block> + Send, R: Reopen<S>, S: SeekableRead {
    bigwig: &'a mut BigWigRead<R, S>,
    known_offset: u64,
    blocks: I,
    vals: Option<Box<dyn Iterator<Item=Value> + Send + 'a>>,
    start: u32,
    end: u32,
}

impl<'a, I, R, S> Iterator for IntervalIter<'a, I, R, S> where I: Iterator<Item=Block> + Send, R: Reopen<S>, S: SeekableRead {
    type Item = io::Result<Value>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match &mut self.vals {
                Some(vals) => {
                    match vals.next() {
                        Some(v) => { return Some(Ok(v)); }
                        None => { self.vals = None; }
                    }
                },
                None => {
                    // TODO: Could minimize this by chunking block reads
                    let current_block = self.blocks.next()?;
                    match get_block_values(self.bigwig, current_block, &mut self.known_offset, self.start, self.end) {
                        Ok(vals) => { self.vals = Some(vals); }
                        Err(e) => { return Some(Err(e)); }
                    }
                },
            }
        }
    }
}

/// Same as IntervalIter but owned
struct OwnedIntervalIter<I, R, S> where I: Iterator<Item=Block> + Send, R: Reopen<S>, S: SeekableRead {
    bigwig: BigWigRead<R, S>,
    known_offset: u64,
    blocks: I,
    vals: Option<Box<dyn Iterator<Item=Value> + Send>>,
    start: u32,
    end: u32,
}

impl<I, R, S> Iterator for OwnedIntervalIter<I, R, S> where I: Iterator<Item=Block> + Send, R: Reopen<S>, S: SeekableRead {
    type Item = io::Result<Value>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match &mut self.vals {
                Some(vals) => {
                    match vals.next() {
                        Some(v) => { return Some(Ok(v)); }
                        None => { self.vals = None; }
                    }
                },
                None => {
                    // TODO: Could minimize this by chunking block reads
                    let current_block = self.blocks.next()?;
                    match get_block_values(&mut self.bigwig, current_block, &mut self.known_offset, self.start, self.end) {
                        Ok(vals) => { self.vals = Some(vals); }
                        Err(e) => { return Some(Err(e)); }
                    }
                },
            }
        }
    }
}

#[derive(Debug)]
pub enum BigWigReadAttachError {
    NotABigWig,
    InvalidChroms,
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

pub struct BigWigRead<R, S> where R:Reopen<S>, S: SeekableRead {
    pub info: BBIFileInfo,
    reopen: R,
    reader: Option<ByteOrdered<BufReader<S>, Endianness>>,
}


impl<R, S> Clone for BigWigRead<R, S> where R: Reopen<S>, S: SeekableRead {
    fn clone(&self) -> Self {
        BigWigRead {
            info: self.info.clone(),
            reopen: self.reopen.clone(),
            reader: None,
        }
    }
}


impl<R: Reopen<S>, S: SeekableRead> BBIRead<S> for BigWigRead<R, S> {
    fn get_info(&self) -> &BBIFileInfo {
        &self.info
    }

    fn ensure_reader(&mut self) -> io::Result<&mut ByteOrdered<BufReader<S>, Endianness>> {
        if self.reader.is_none() {
            let endianness = self.info.header.endianness;
            let fp = self.reopen.reopen()?;
            let file = ByteOrdered::runtime(BufReader::new(fp), endianness);
            self.reader.replace(file);
        }
        Ok(self.reader.as_mut().unwrap())
    }

    fn close(&mut self) {
        self.reader.take();
    }

    fn get_chroms(&self) -> Vec<ChromAndSize> {
        self.info.chrom_info.iter().map(|c| ChromAndSize { name: c.name.clone(), length: c.length }).collect::<Vec<_>>()
    }
}

impl BigWigRead<ReopenableFile, File> {
    pub fn from_file_and_attach(path: String) -> Result<Self, BigWigReadAttachError> {
        let reopen = ReopenableFile { path: path.clone() };
        let b = BigWigRead::from(reopen);
        if b.is_err() {
            eprintln!("Error when opening: {}", path);
        }
        b
    }
}

impl<R,S> BigWigRead<R,S> where R: Reopen<S>, S: SeekableRead {
    pub fn from(reopen: R) -> Result<Self, BigWigReadAttachError> {
        let fp = reopen.reopen()?;
        let file = BufReader::new(fp);
        let info = match read_info(file) {
            Err(e) => {
                return Err(e.into());
            }
            Ok(info) => info,
        };
        match info.filetype {
            BBIFile::BigWig => {},
            _ => return Err(BigWigReadAttachError::NotABigWig),
        }

        Ok(BigWigRead {
            info,
            reopen,
            reader: None,
        })
    }
}

impl<R: Reopen<S> + 'static, S: SeekableRead + 'static> BigWigRead<R, S> {
    pub fn get_interval_move(mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<impl Iterator<Item=io::Result<Value>> + Send> {
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
        Ok(OwnedIntervalIter {
            bigwig: self,
            known_offset: 0,
            blocks: blocks.into_iter(),
            vals: None,
            start,
            end,
        })
    }

    pub fn get_interval<'a>(&'a mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<impl Iterator<Item=io::Result<Value>> + Send + 'a> {
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
        Ok(IntervalIter {
            bigwig: self,
            known_offset: 0,
            blocks: blocks.into_iter(),
            vals: None,
            start,
            end,
        })
    }

    pub fn get_zoom_interval<'a>(&'a mut self, chrom_name: &str, start: u32, end: u32, reduction_level: u32) -> io::Result<impl Iterator<Item=io::Result<ZoomRecord>> + Send + 'a> {
        let zoom_header = match self.info.zoom_headers.iter().find(|h| h.reduction_level == reduction_level) {
            Some(h) => h,
            None => {
                return Err(io::Error::new(io::ErrorKind::Other, "No reduction level found."));
            }
        };

        let index_offset = zoom_header.index_offset;
        let file = self.ensure_reader()?;
        file.seek(SeekFrom::Start(index_offset))?;
        let blocks = self.search_cir_tree(chrom_name, start, end)?;
        Ok(ZoomIntervalIter::new(self, blocks.into_iter(), start, end))
    }
}

fn get_block_values<R: Reopen<S>, S: SeekableRead>(bigwig: &mut BigWigRead<R, S>, block: Block, known_offset: &mut u64, start: u32, end: u32) -> io::Result<Box<dyn Iterator<Item=Value> + Send>> {
    let mut block_data_mut = get_block_data(bigwig, &block, *known_offset)?;
    let mut values: Vec<Value> = Vec::new();

    let _chrom_id = block_data_mut.read_u32()?;
    let chrom_start = block_data_mut.read_u32()?;
    let _chrom_end = block_data_mut.read_u32()?;
    let item_step = block_data_mut.read_u32()?;
    let item_span = block_data_mut.read_u32()?;
    let section_type = block_data_mut.read_u8()?;
    let _reserved = block_data_mut.read_u8()?;
    let item_count = block_data_mut.read_u16()?;

    let mut curr_start = chrom_start;
    for _ in 0..item_count {
        let mut value = match section_type {
            1 => {
                // bedgraph
                let chrom_start = block_data_mut.read_u32()?;
                let chrom_end = block_data_mut.read_u32()?;
                let value = block_data_mut.read_f32()?;
                Value {
                    start: chrom_start,
                    end: chrom_end,
                    value,
                }
            },
            2 => {
                // variable step
                let chrom_start = block_data_mut.read_u32()?;
                let chrom_end = chrom_start + item_span;
                let value = block_data_mut.read_f32()?;
                Value {
                    start: chrom_start,
                    end: chrom_end,
                    value,
                }
            },
            3 => {
                // fixed step
                let chrom_start = curr_start;
                curr_start += item_step;
                let chrom_end = chrom_start + item_span;
                let value = block_data_mut.read_f32()?;
                Value {
                    start: chrom_start,
                    end: chrom_end,
                    value,
                }
            },
            _ => return Err(std::io::Error::new(std::io::ErrorKind::Other, format!("Unknown bigwig section type: {}", section_type)))
        };
        if value.end >= start && value.start <= end {
            value.start = value.start.max(start);
            value.end = value.end.min(end);
            values.push(value)
        }
    }

    *known_offset = block.offset + block.size;
    Ok(Box::new(values.into_iter()))
}
