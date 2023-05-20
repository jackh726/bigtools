use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::vec::Vec;

use byteordered::{ByteOrdered, Endianness};
use thiserror::Error;

use crate::bbiread::{
    get_block_data, read_info, BBIFileInfo, BBIFileReadInfoError, BBIRead, BBIReadError, BBIReader,
    Block, ChromAndSize, MemCachedReader, ZoomIntervalIter,
};
use crate::bigwig::{BBIFile, BedEntry, ZoomRecord};
use crate::utils::mem_cached_file::{MemCachedRead, CACHE_SIZE};
use crate::utils::seekableread::{Reopen, ReopenableFile, SeekableRead};

struct IntervalIter<'a, I, R, S>
where
    I: Iterator<Item = Block> + Send,
    R: Reopen<S>,
    S: SeekableRead,
{
    bigbed: &'a mut BigBedRead<R, S>,
    known_offset: u64,
    blocks: I,
    vals: Option<Box<dyn Iterator<Item = BedEntry> + Send + 'a>>,
    expected_chrom: u32,
    start: u32,
    end: u32,
}

impl<'a, I, R, S> Iterator for IntervalIter<'a, I, R, S>
where
    I: Iterator<Item = Block> + Send,
    R: Reopen<S>,
    S: SeekableRead,
{
    type Item = Result<BedEntry, BBIReadError>;

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
                    match get_block_entries(
                        self.bigbed,
                        current_block,
                        &mut self.known_offset,
                        self.expected_chrom,
                        self.start,
                        self.end,
                    ) {
                        Ok(vals) => {
                            self.vals = Some(vals);
                        }
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
    bigbed: BigBedRead<R, S>,
    known_offset: u64,
    blocks: I,
    vals: Option<Box<dyn Iterator<Item = BedEntry> + Send>>,
    expected_chrom: u32,
    start: u32,
    end: u32,
}

impl<I, R, S> Iterator for OwnedIntervalIter<I, R, S>
where
    I: Iterator<Item = Block> + Send,
    R: Reopen<S>,
    S: SeekableRead,
{
    type Item = Result<BedEntry, BBIReadError>;

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
                    match get_block_entries(
                        &mut self.bigbed,
                        current_block,
                        &mut self.known_offset,
                        self.expected_chrom,
                        self.start,
                        self.end,
                    ) {
                        Ok(vals) => {
                            self.vals = Some(vals);
                        }
                        Err(e) => {
                            return Some(Err(e));
                        }
                    }
                }
            }
        }
    }
}

#[derive(Error, Debug)]
pub enum BigBedReadAttachError {
    #[error("File is not a bigBed.")]
    NotABigBed,
    #[error("The chromosomes are invalid.")]
    InvalidChroms,
    #[error("An error occurred: {}", .0)]
    IoError(io::Error),
}

impl From<io::Error> for BigBedReadAttachError {
    fn from(error: io::Error) -> Self {
        BigBedReadAttachError::IoError(error)
    }
}

impl From<BBIFileReadInfoError> for BigBedReadAttachError {
    fn from(error: BBIFileReadInfoError) -> Self {
        match error {
            BBIFileReadInfoError::UnknownMagic => BigBedReadAttachError::NotABigBed,
            BBIFileReadInfoError::InvalidChroms => BigBedReadAttachError::InvalidChroms,
            BBIFileReadInfoError::IoError(e) => BigBedReadAttachError::IoError(e),
        }
    }
}

pub struct BigBedRead<R, S> {
    pub info: BBIFileInfo,
    reopen: R,
    reader: Option<ByteOrdered<BufReader<S>, Endianness>>,
    cache: HashMap<usize, [u8; CACHE_SIZE]>,
}

impl<R, S> Clone for BigBedRead<R, S>
where
    R: Reopen<S>,
    S: SeekableRead,
{
    fn clone(&self) -> Self {
        BigBedRead {
            info: self.info.clone(),
            reopen: self.reopen.clone(),
            reader: None,
            cache: HashMap::new(),
        }
    }
}

impl<R, S> BBIRead<S> for BigBedRead<R, S>
where
    R: Reopen<S>,
    S: SeekableRead,
{
    fn get_info(&self) -> &BBIFileInfo {
        &self.info
    }

    fn autosql(&mut self) -> io::Result<String> {
        self.ensure_reader()?;
        let reader = self.reader.as_mut().unwrap();
        reader.seek(SeekFrom::Start(self.info.header.auto_sql_offset))?;
        let mut buffer = Vec::new();
        reader.read_until(b'\0', &mut buffer)?;
        buffer.pop();
        let autosql = String::from_utf8(buffer)
            .map_err(|_| io::Error::new(io::ErrorKind::Other, "Invalid autosql: not UTF-8"))?;
        Ok(autosql)
    }

    fn ensure_reader(&mut self) -> io::Result<&mut BBIReader<S>> {
        if self.reader.is_none() {
            let endianness = self.info.header.endianness;
            let fp = self.reopen.reopen()?;
            let file = ByteOrdered::runtime(BufReader::new(fp), endianness);
            self.reader.replace(file);
        }
        // FIXME: In theory, can get rid of this unwrap by doing a `match` with
        // `Option::insert` in the `None` case, but that currently runs into
        // lifetime issues.
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

impl BigBedRead<ReopenableFile, File> {
    pub fn from_file_and_attach(path: String) -> Result<Self, BigBedReadAttachError> {
        let reopen = ReopenableFile { path: path.clone() };
        let b = BigBedRead::from(reopen);
        if b.is_err() {
            eprintln!("Error when opening: {}", path);
        }
        b
    }
}

impl<R, S> BigBedRead<R, S>
where
    R: Reopen<S>,
    S: SeekableRead,
{
    pub fn from(reopen: R) -> Result<Self, BigBedReadAttachError> {
        let file = reopen.reopen()?;
        let info = read_info(file)?;
        match info.filetype {
            BBIFile::BigBed => {}
            _ => return Err(BigBedReadAttachError::NotABigBed),
        }

        Ok(BigBedRead {
            info,
            reopen,
            reader: None,
            cache: HashMap::new(),
        })
    }

    pub fn get_interval<'a>(
        &'a mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> Result<impl Iterator<Item = Result<BedEntry, BBIReadError>> + Send + 'a, BBIReadError>
    {
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
        // TODO: this is only for asserting that the chrom is what we expect
        let chrom_ix = self
            .get_info()
            .chrom_info
            .iter()
            .find(|&x| x.name == chrom_name)
            .unwrap()
            .id;
        Ok(IntervalIter {
            bigbed: self,
            known_offset: 0,
            blocks: blocks.into_iter(),
            vals: None,
            expected_chrom: chrom_ix,
            start,
            end,
        })
    }

    pub fn get_interval_move(
        mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> Result<impl Iterator<Item = Result<BedEntry, BBIReadError>> + Send, BBIReadError> {
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
        // TODO: this is only for asserting that the chrom is what we expect
        let chrom_ix = self
            .get_info()
            .chrom_info
            .iter()
            .find(|&x| x.name == chrom_name)
            .unwrap()
            .id;
        Ok(OwnedIntervalIter {
            bigbed: self,
            known_offset: 0,
            blocks: blocks.into_iter(),
            vals: None,
            expected_chrom: chrom_ix,
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
}

// TODO: remove expected_chrom
fn get_block_entries<R: Reopen<S>, S: SeekableRead>(
    bigbed: &mut BigBedRead<R, S>,
    block: Block,
    known_offset: &mut u64,
    expected_chrom: u32,
    start: u32,
    end: u32,
) -> Result<Box<dyn Iterator<Item = BedEntry> + Send>, BBIReadError> {
    let mut block_data_mut = get_block_data(bigbed, &block, *known_offset)?;
    let mut entries: Vec<BedEntry> = Vec::new();

    let mut read_entry = || -> io::Result<BedEntry> {
        let chrom_id = block_data_mut.read_u32()?;
        let chrom_start = block_data_mut.read_u32()?;
        let chrom_end = block_data_mut.read_u32()?;
        if chrom_start == 0 && chrom_end == 0 {
            return Err(io::Error::new(io::ErrorKind::Other, ""));
        }
        // FIXME: should this just return empty?
        assert_eq!(
            chrom_id, expected_chrom,
            "BUG: bigBed had multiple chroms in a section"
        );
        let s: Vec<u8> = block_data_mut
            .by_ref()
            .bytes()
            .take_while(|c| {
                if let Ok(c) = c {
                    return *c != b'\0';
                }
                false
            })
            .collect::<Result<Vec<u8>, _>>()?;
        let rest = String::from_utf8(s).unwrap();
        Ok(BedEntry {
            start: chrom_start,
            end: chrom_end,
            rest,
        })
    };
    while let Ok(entry) = read_entry() {
        if entry.end >= start && entry.start <= end {
            entries.push(entry);
        }
    }

    *known_offset = block.offset + block.size;
    Ok(Box::new(entries.into_iter()))
}
