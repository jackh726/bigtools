use std::io::{self, Read, Seek, SeekFrom};
use std::io::{BufReader};
use std::fs::File;
use std::vec::Vec;

use byteordered::{ByteOrdered, Endianness};
use flate2::read::ZlibDecoder;

use crate::seekableread::SeekableRead;
use crate::bbiread::{BBIRead, BBIFileReadInfoError, BBIFileInfo, Block, ChromAndSize};
use crate::bigwig::{BBIFile, BedEntry};


fn get_entries<R: SeekableRead + 'static>(bigbed: &mut BigBedRead<R>, current_block: Block, known_offset: &mut u64, expected_chrom: u32) -> io::Result<Box<dyn Iterator<Item=BedEntry> + Send>> {
    let endianness = bigbed.info.header.endianness;
    let uncompress_buf_size: usize = bigbed.info.header.uncompress_buf_size as usize;
    let file = bigbed.reader.as_mut().unwrap();

    // TODO: Could minimize this by chunking block reads
    if *known_offset != current_block.offset {
        match file.seek(SeekFrom::Start(current_block.offset)) {
            Ok(_) => {},
            Err(e) => return Err(e),
        }
    }
    let block_vals = match BigBedRead::get_block_entries(file, &current_block, endianness, uncompress_buf_size, expected_chrom) {
        Ok(vals) => vals,
        Err(e) => return Err(e),
    };
    *known_offset = current_block.offset + current_block.size;
    Ok(Box::new(block_vals))
}

struct IntervalIter<'a, I, R: SeekableRead> where I: Iterator<Item=Block> + Send {
    bigbed: &'a mut BigBedRead<R>,
    known_offset: u64,
    blocks: I,
    vals: Option<Box<dyn Iterator<Item=BedEntry> + Send + 'a>>,
    expected_chrom: u32,
}

impl<'a, I, R: SeekableRead + 'static> Iterator for IntervalIter<'a, I, R> where I: Iterator<Item=Block> + Send {
    type Item = io::Result<BedEntry>;

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
                    match get_entries(self.bigbed, current_block, &mut self.known_offset, self.expected_chrom) {
                        Ok(vals) => { self.vals = Some(vals); }
                        Err(e) => { return Some(Err(e)); }
                    }
                },
            }
        }
    }
}

/// Same as IntervalIter but owned
struct OwnedIntervalIter<I, R: SeekableRead> where I: Iterator<Item=Block> + Send {
    bigbed: BigBedRead<R>,
    known_offset: u64,
    blocks: I,
    vals: Option<Box<dyn Iterator<Item=BedEntry> + Send>>,
    expected_chrom: u32,
}

impl<I, R: SeekableRead + 'static> Iterator for OwnedIntervalIter<I, R> where I: Iterator<Item=Block> + Send {
    type Item = io::Result<BedEntry>;

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
                    match get_entries(&mut self.bigbed, current_block, &mut self.known_offset, self.expected_chrom) {
                        Ok(vals) => { self.vals = Some(vals); }
                        Err(e) => { return Some(Err(e)); }
                    }
                },
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
        return match error {
            BBIFileReadInfoError::UnknownMagic => BigBedReadAttachError::NotABigBed,
            BBIFileReadInfoError::InvalidChroms => BigBedReadAttachError::InvalidChroms,
            BBIFileReadInfoError::IoError(e) => BigBedReadAttachError::IoError(e),
        }
    }
}

pub struct BigBedRead<R: SeekableRead> {
    pub path: String,
    pub info: BBIFileInfo,
    reader: Option<ByteOrdered<BufReader<R>, Endianness>>,
}

impl<R: SeekableRead> Clone for BigBedRead<R> {
    fn clone(&self) -> Self {
        BigBedRead {
            path: self.path.clone(),
            info: self.info.clone(),
            reader: None,
        }
    }
}

impl BBIRead<File> for BigBedRead<File> {
    fn get_info(&self) -> &BBIFileInfo {
        &self.info
    }

    fn ensure_reader(&mut self) -> io::Result<&mut ByteOrdered<BufReader<File>, Endianness>> {
        if self.reader.is_none() {
            let endianness = self.info.header.endianness;
            let fp = File::open(self.path.clone())?;
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

impl BigBedRead<File> {
    pub fn from_file_and_attach(path: String) -> Result<Self, BigBedReadAttachError> {
        let fp = File::open(path.clone())?;
        let file = BufReader::new(fp);
        let info = match BigBedRead::read_info(file) {
            Err(e) => {
                eprintln!("Error when opening: {}", path.clone());
                return Err(e.into());
            }
            Ok(info) => info,
        };
        match info.filetype {
            BBIFile::BigBed => {},
            _ => return Err(BigBedReadAttachError::NotABigBed),
        }

        Ok(BigBedRead {
            path,
            info,
            reader: None,
        })
    }

    // TODO: remove expected_chrom
    /// This assumes that the file is currently at the block's start
    fn get_block_entries<R: SeekableRead>(file: &mut ByteOrdered<BufReader<R>, Endianness>, block: &Block, endianness: Endianness, uncompress_buf_size: usize, expected_chrom: u32) -> io::Result<impl Iterator<Item=BedEntry>> {
        let mut entries: Vec<BedEntry> = Vec::new();

        let mut raw_data = vec![0u8; block.size as usize];
        file.read_exact(&mut raw_data)?;
        let block_data: Vec<u8> = if uncompress_buf_size > 0 {
            let mut uncompressed_block_data = vec![0u8; uncompress_buf_size];
            let mut d = ZlibDecoder::new(&raw_data[..]);
            let _ = d.read(&mut uncompressed_block_data)?;
            uncompressed_block_data
        } else {
            raw_data
        };

        let mut block_data_mut = ByteOrdered::runtime(&block_data[..], endianness);

        let mut read_entry = || -> io::Result<BedEntry> {
            let _chrom_id = block_data_mut.read_u32()?;
            assert_eq!(_chrom_id, expected_chrom, "BUG: bigBed had multiple chroms in a section");
            let chrom_start = block_data_mut.read_u32()?;
            let chrom_end = block_data_mut.read_u32()?;
            let s: Vec<u8> = block_data_mut.by_ref().bytes().take_while(|c| {
                if let Ok(c) = c {
                    return *c != b'\0';
                }
                false
            }).collect::<Result<Vec<u8>,_>>()?;
            let rest = String::from_utf8(s).unwrap();
            Ok(BedEntry {
                start: chrom_start,
                end: chrom_end,
                rest: rest,
            })
        };
        loop {
            match read_entry() {
                Ok(entry) => {
                    // TODO: the entire section could be terminated by many 0s. Need to identify a better way of filtering out these    
                    if entry.start == 0 && entry.end == 0 {
                        break
                    }
                    entries.push(entry)
                },
                Err(_) => break,
            }
        }

        Ok(entries.into_iter())
    }

    pub fn get_interval<'a>(&'a mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<impl Iterator<Item=io::Result<BedEntry>> + Send + 'a> {
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
        // TODO: this is only for asserting that the chrom is what we expect
        let chrom_ix = self.get_info().chrom_info.iter().find(|&x| x.name == chrom_name).unwrap().id;
        Ok(IntervalIter {
            bigbed: self,
            known_offset: 0,
            blocks: blocks.into_iter(),
            vals: None,
            expected_chrom: chrom_ix,
        })
    }

    pub fn get_interval_move(mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<impl Iterator<Item=io::Result<BedEntry>> + Send> {
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
        // TODO: this is only for asserting that the chrom is what we expect
        let chrom_ix = self.get_info().chrom_info.iter().find(|&x| x.name == chrom_name).unwrap().id;
        Ok(OwnedIntervalIter {
            bigbed: self,
            known_offset: 0,
            blocks: blocks.into_iter(),
            vals: None,
            expected_chrom: chrom_ix,
        })
    }

}