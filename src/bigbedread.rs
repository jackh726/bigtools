use std::io::{self, Read, Seek, SeekFrom};
use std::io::{BufReader};
use std::fs::File;
use std::vec::Vec;

use byteordered::{ByteOrdered, Endianness};
use flate2::read::ZlibDecoder;

use crate::bbiread::{BBIFile, BBIRead, BBIFileReadInfoError, BBIFileInfo, Block, ChromAndSize};
use crate::bigwig::BedEntry;


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

pub struct BigBedRead {
    pub path: String,
    pub info: BBIFileInfo,
    reader: Option<ByteOrdered<BufReader<File>, Endianness>>,
}

impl Clone for BigBedRead {
    fn clone(&self) -> Self {
        BigBedRead {
            path: self.path.clone(),
            info: self.info.clone(),
            reader: None,
        }
    }
}

impl BBIRead for BigBedRead {
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

impl BigBedRead {
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

    /// This assumes that the file is currently at the block's start
    fn get_block_entries(file: &mut ByteOrdered<BufReader<File>, Endianness>, block: &Block, endianness: Endianness, uncompress_buf_size: usize, expected_chrom: u32) -> io::Result<impl Iterator<Item=BedEntry>> {
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

        let (endianness, uncompress_buf_size) = {
            let info = self.get_info();
            (info.header.endianness, info.header.uncompress_buf_size as usize)
        };

        let mut file = self.ensure_reader()?;

        if blocks.len() > 0 {
            file.seek(SeekFrom::Start(blocks[0].offset))?;
        }
        let mut iter = blocks.into_iter().peekable();
        
        let block_iter = std::iter::from_fn(move || {
            let next = iter.next();
            let peek = iter.peek();
            next.map(|n| (n, peek.map(|p| p.offset)))
        });
        let vals_iter = block_iter.flat_map(move |(block, next_offset)| {
            let mut vals = || -> io::Result<Box<dyn Iterator<Item=io::Result<BedEntry>> + Send + 'a>> {
                // TODO: Could minimize this by chunking block reads
                let vals = BigBedRead::get_block_entries(&mut file, &block, endianness, uncompress_buf_size, chrom_ix)?;
                if let Some(next_offset) = next_offset {
                    if next_offset != block.offset + block.size {
                        file.seek(SeekFrom::Start(next_offset))?;
                    }
                }
                Ok(Box::new(vals.map(|v| Ok(v))))
            };
            let v: Box<dyn Iterator<Item=io::Result<BedEntry>> + Send + 'a> = vals().unwrap_or_else(|e| Box::new(std::iter::once(Err(e))));
            v
        }).filter_map(move |mut val| {
            if let Ok(ref mut v) = val {
                if v.end < start || v.start >= end {
                    return None;
                }
            }
            return Some(val);
        });

        Ok(vals_iter)
    }

}