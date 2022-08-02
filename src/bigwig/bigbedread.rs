use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::vec::Vec;

use byteordered::ByteOrdered;

use crate::bbiread::{
    get_block_data, read_info, BBIFileInfo, BBIFileReadInfoError, BBIRead, BBIReader, Block,
    ChromAndSize, MemCachedReader, ZoomIntervalIter,
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
    type Item = io::Result<BedEntry>;

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
    type Item = io::Result<BedEntry>;

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

#[derive(Debug)]
pub enum BigBedReadAttachError {
    NotABigBed,
    InvalidChroms,
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
    reader: Option<S>,
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
        let auto_sql_offset = self.info.header.auto_sql_offset;
        let mut reader = self.ensure_reader()?;
        reader.seek(SeekFrom::Start(auto_sql_offset))?;
        let mut buffer = Vec::new();
        reader.read_until(b'\0', &mut buffer)?;
        buffer.pop();
        let autosql = String::from_utf8(buffer)
            .map_err(|_| io::Error::new(io::ErrorKind::Other, "Invalid autosql: not UTF-8"))?;
        Ok(autosql)
    }

    fn ensure_reader(&mut self) -> io::Result<BBIReader<&mut S>> {
        if self.reader.is_none() {
            let fp = self.reopen.reopen()?;
            self.reader.replace(fp);
        }
        // FIXME: In theory, can get rid of this unwrap by doing a `match` with
        // `Option::insert` in the `None` case, but that currently runs into
        // lifetime issues.
        let endianness = self.info.header.endianness;
        let reader =
            ByteOrdered::runtime(BufReader::new(self.reader.as_mut().unwrap()), endianness);
        Ok(reader)
    }

    fn ensure_mem_cached_reader(&mut self) -> io::Result<MemCachedReader<'_, S>> {
        if self.reader.is_none() {
            let fp = self.reopen.reopen()?;
            self.reader.replace(fp);
        }
        let endianness = self.info.header.endianness;
        let inner = self.reader.as_mut().unwrap();
        Ok(ByteOrdered::runtime(
            MemCachedRead::new(inner, &mut self.cache),
            endianness,
        ))
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
    ) -> io::Result<impl Iterator<Item = io::Result<BedEntry>> + Send + 'a> {
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
    ) -> io::Result<impl Iterator<Item = io::Result<BedEntry>> + Send> {
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
    ) -> io::Result<impl Iterator<Item = io::Result<ZoomRecord>> + Send + 'a> {
        let chrom = self.info.chrom_id(chrom_name)?;
        let zoom_header = match self
            .info
            .zoom_headers
            .iter()
            .find(|h| h.reduction_level == reduction_level)
        {
            Some(h) => h,
            None => {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "No reduction level found.",
                ));
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
    expected_chrom: u32,
    start: u32,
    end: u32,
) -> io::Result<Box<dyn Iterator<Item = BedEntry> + Send>> {
    let mut block_data_mut = get_block_data(bigbed, &block)?;
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

    Ok(Box::new(entries.into_iter()))
}
