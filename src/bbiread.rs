use std::io::{self, Read, Seek, SeekFrom};
use std::io::{BufReader};
use std::fs::File;
use std::vec::Vec;

use byteordered::{ByteOrdered, Endianness};

use crate::bigwig::{BBIFile, ZoomHeader, CHROM_TREE_MAGIC, CIR_TREE_MAGIC, BIGWIG_MAGIC_LTH, BIGWIG_MAGIC_HTL, BIGBED_MAGIC_LTH, BIGBED_MAGIC_HTL};


#[derive(Debug)]
pub struct Block {
    pub offset: u64,
    pub size: u64,
}

#[derive(Clone, Debug)]
pub struct BBIHeader {
    pub endianness: Endianness,

    pub(crate) version: u16,
    pub(crate) zoom_levels: u16,
    pub(crate) chromosome_tree_offset: u64,
    pub(crate) full_data_offset: u64,
    pub(crate) full_index_offset: u64,
    pub(crate) field_count: u16,
    pub(crate) defined_field_count: u16,
    pub(crate) auto_sql_offset: u64,
    pub(crate) total_summary_offset: u64,
    pub(crate) uncompress_buf_size: u32,
}

#[derive(Clone, Debug)]
pub struct ChromInfo {
    pub name: String,
    pub length: u32,
    pub(crate) id: u32,
}

#[derive(Debug)]
pub struct ChromAndSize {
    pub name: String,
    pub length: u32,
}

impl PartialEq for ChromAndSize {
    fn eq(&self, other: &ChromAndSize) -> bool {
        self.name == other.name
    }
}

#[derive(Clone, Debug)]
pub struct BBIFileInfo {
    pub filetype: BBIFile,
    pub header: BBIHeader,
    pub(crate) zoom_headers: Vec<ZoomHeader>,
    pub(crate) chrom_info: Vec<ChromInfo>,
}

pub enum BBIFileReadInfoError {
    UnknownMagic,
    InvalidChroms,
    IoError(io::Error),
}

impl From<io::Error> for BBIFileReadInfoError {
    fn from(error: io::Error) -> Self {
        BBIFileReadInfoError::IoError(error)
    }
}

pub trait BBIRead {
    fn get_info(&self) -> &BBIFileInfo;

    fn ensure_reader(&mut self) -> io::Result<&mut ByteOrdered<BufReader<File>, Endianness>>;

    /// Manually close the open file descriptor (if it exists). If any operations are performed after this is called, the file descriptor will be reopened.
    fn close(&mut self);

    fn get_chroms(&self) -> Vec<ChromAndSize>;


    fn read_info(file: BufReader<File>) -> Result<BBIFileInfo, BBIFileReadInfoError> {
        let mut file = ByteOrdered::runtime(file, Endianness::Little);

        let magic = file.read_u32()?;
        // println!("Magic {:x?}", magic);
        let filetype = match magic {
            BIGWIG_MAGIC_HTL => {
                file = file.into_opposite();
                BBIFile::BigWig
            },
            BIGWIG_MAGIC_LTH => {
                BBIFile::BigWig
            },
            BIGBED_MAGIC_HTL => {
                file = file.into_opposite();
                BBIFile::BigBed
            },
            BIGBED_MAGIC_LTH => {
                BBIFile::BigBed
            },
            _ => return Err(BBIFileReadInfoError::UnknownMagic),
        };

        let version = file.read_u16()?;

        // TODO: should probably handle versions < 3
        let zoom_levels = file.read_u16()?;
        let chromosome_tree_offset = file.read_u64()?;
        let full_data_offset = file.read_u64()?;
        let full_index_offset = file.read_u64()?;
        let field_count = file.read_u16()?;
        let defined_field_count = file.read_u16()?;
        let auto_sql_offset = file.read_u64()?;
        let total_summary_offset = file.read_u64()?;
        let uncompress_buf_size = file.read_u32()?;
        let _reserved = file.read_u64()?;

        let header = BBIHeader {
            endianness: file.endianness(),
            version,
            zoom_levels,
            chromosome_tree_offset,
            full_data_offset,
            full_index_offset,
            field_count,
            defined_field_count,
            auto_sql_offset,
            total_summary_offset,
            uncompress_buf_size,
        };

        //println!("Header: {:?}", header);

        let zoom_headers = read_zoom_headers(&mut file, &header)?;

        // TODO: could instead store this as an Option and only read when needed
        file.seek(SeekFrom::Start(header.chromosome_tree_offset))?;
        let magic = file.read_u32()?;
        let _block_size = file.read_u32()?;
        let key_size = file.read_u32()?;
        let val_size = file.read_u32()?;
        let item_count = file.read_u64()?;
        let _reserved = file.read_u64()?;
        if magic != CHROM_TREE_MAGIC {
            return Err(BBIFileReadInfoError::InvalidChroms);
        }
        //println!("{:x?} {:?} {:?} {:?} {:?} {:?}", magic, _block_size, key_size, val_size, item_count, _reserved);
        assert_eq!(val_size, 8u32); 

        let mut chrom_info = Vec::with_capacity(item_count as usize);
        read_chrom_tree_block(&mut file, &mut chrom_info, key_size)?;

        let info = BBIFileInfo {
            filetype,
            header,
            zoom_headers,
            chrom_info,
        };

        Ok(info)
    }

    /// This assumes the file is at the cir tree start
    fn search_cir_tree(&mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<Vec<Block>> {
        let chrom_ix = {
            let chrom_info = &self.get_info().chrom_info;
            let chrom = chrom_info.iter().find(|&x| x.name == chrom_name);
            //println!("Chrom: {:?}", chrom);
            match chrom {
                Some(c) => c.id,
                None => return Err(std::io::Error::new(std::io::ErrorKind::Other, format!("{} not found.", chrom_name)))
            }
        };

        let mut file = self.ensure_reader()?;

        let magic = file.read_u32()?;
        if magic != CIR_TREE_MAGIC {
            return Err(std::io::Error::new(std::io::ErrorKind::Other, "Invalid file format: CIR_TREE_MAGIC does not match."));
        }
        let _blocksize = file.read_u32()?;
        let _item_count = file.read_u64()?;
        let _start_chrom_idx = file.read_u32()?;
        let _start_base = file.read_u32()?;
        let _end_chrom_idx = file.read_u32()?;
        let _end_base = file.read_u32()?;
        let _end_file_offset = file.read_u64()?;
        let _item_per_slot = file.read_u32()?;
        let _reserved = file.read_u32()?;

        // TODO: could do some optimization here to check if our interval overlaps with any data

        //println!("cirTree header:\n bs: {:?}\n ic: {:?}\n sci: {:?}\n sb: {:?}\n eci: {:?}\n eb: {:?}\n efo: {:?}\n ips: {:?}\n r: {:?}", _blocksize, _item_count, _start_chrom_idx, _start_base, _end_chrom_idx, _end_base, _end_file_offset, _item_per_slot, _reserved);
        let mut blocks: Vec<Block> = vec![];
        search_overlapping_blocks(&mut file, chrom_ix, start, end, &mut blocks)?;
        //println!("overlapping_blocks: {:?}", blocks);
        Ok(blocks)
    }

    fn get_overlapping_blocks(&mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<Vec<Block>> {
        let full_index_offset = self.get_info().header.full_index_offset;

        let file = self.ensure_reader()?;
        file.seek(SeekFrom::Start(full_index_offset))?;

        self.search_cir_tree(chrom_name, start, end)
    }
}

fn read_zoom_headers(file: &mut ByteOrdered<BufReader<File>, Endianness>, header: &BBIHeader) -> io::Result<Vec<ZoomHeader>> {
    let mut zoom_headers = vec![];
    for _ in 0..header.zoom_levels {
        let reduction_level = file.read_u32()?;
        let _reserved = file.read_u32()?;
        let data_offset = file.read_u64()?;
        let index_offset = file.read_u64()?;

        //println!("Zoom header: reductionLevel: {:?} Reserved: {:?} Data offset: {:?} Index offset: {:?}", reduction_level, _reserved, data_offset, index_offset);

        zoom_headers.push(ZoomHeader {
            reduction_level,
            data_offset,
            index_offset,
        });
    }

    Ok(zoom_headers)
}

fn read_chrom_tree_block(f: &mut ByteOrdered<BufReader<File>, Endianness>, chroms: &mut Vec<ChromInfo>, key_size: u32) -> io::Result<()> {
    let isleaf = f.read_u8()?;
    let _reserved = f.read_u8()?;
    let count = f.read_u16()?;

    if isleaf == 1 {
        for _ in 0..count {
            let mut key_bytes = vec![0u8; key_size as usize];
            f.read_exact(&mut key_bytes)?;
            let key_string = match String::from_utf8(key_bytes) {
                Ok(s) => s.trim_matches(char::from(0)).to_owned(),
                Err(_) => return Err(io::Error::new(io::ErrorKind::Other, "Invalid file format: Invalid utf-8 string.")),
            };
            let chrom_id = f.read_u32()?;
            let chrom_size = f.read_u32()?;
            chroms.push(ChromInfo {
                name: key_string,
                id: chrom_id,
                length: chrom_size,
            });
        }
    } else {
        let mut current_position: u64;
        for _ in 0..count {
            let mut key_bytes = vec![0u8; key_size as usize];
            f.read_exact(&mut key_bytes)?;
            // TODO: could add specific find here by comparing key string
            let child_offset = f.read_u64()?;
            current_position = f.seek(SeekFrom::Current(0))?;
            f.seek(SeekFrom::Start(child_offset))?;
            read_chrom_tree_block(f, chroms, key_size)?;
            f.seek(SeekFrom::Start(current_position))?;
        }
    }
    Ok(())
}

#[inline]
fn compare_position(chrom1: u32, chrom1_base: u32, chrom2: u32, chrom2_base: u32) -> i8 {
    if chrom1 < chrom2 {
        -1
    } else if chrom1 > chrom2 {
        1
    } else if chrom1_base < chrom2_base {
        -1
    } else if chrom1_base > chrom2_base {
        1
    } else {
        0
    }
}

#[inline]
fn overlaps(chromq: u32, chromq_start: u32, chromq_end: u32, chromb1: u32, chromb1_start: u32, chromb2: u32, chromb2_end: u32) -> bool {
    compare_position(chromq, chromq_start, chromb2, chromb2_end) <= 0 && compare_position(chromq, chromq_end, chromb1, chromb1_start) >= 0
}

fn search_overlapping_blocks(mut file: &mut ByteOrdered<BufReader<File>, Endianness>, chrom_ix: u32, start: u32, end: u32, mut blocks: &mut Vec<Block>) -> io::Result<()> {
    //println!("Searching for overlapping blocks at {:?}. Searching {:?}:{:?}-{:?}", self.current_file_offset()?, chrom_ix, start, end);

    let isleaf: u8 = file.read_u8()?;
    assert!(isleaf == 1 || isleaf == 0, "Unexpected isleaf: {}", isleaf);
    let _reserved = file.read_u8()?;
    let count: u16 = file.read_u16()?;
    //println!("Index: {:?} {:?} {:?}", isleaf, _reserved, count);

    let mut childblocks: Vec<u64> = vec![];
    for _ in 0..count {
        let start_chrom_ix = file.read_u32()?;
        let start_base = file.read_u32()?;
        let end_chrom_ix = file.read_u32()?;
        let end_base = file.read_u32()?;
        if isleaf == 1 {
            let data_offset = file.read_u64()?;
            let data_size = file.read_u64()?;
            if !overlaps(chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base) {
                continue;
            }
            println!("Overlaps (leaf): {:?}:{:?}-{:?} with {:?}:{:?}-{:?}:{:?} {:?} {:?}", chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset, data_size);
            blocks.push(Block {
                offset: data_offset,
                size: data_size,
            })
        } else {
            let data_offset = file.read_u64()?;
            if !overlaps(chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base) {
                continue;
            }
            println!("Overlaps (non-leaf): {:?}:{:?}-{:?} with {:?}:{:?}-{:?}:{:?} {:?}", chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset);
            childblocks.push(data_offset);
        }
    }
    for childblock in childblocks {
        //println!("Seeking to {:?}", childblock);
        file.seek(SeekFrom::Start(childblock))?;
        search_overlapping_blocks(&mut file, chrom_ix, start, end, &mut blocks)?;
    }
    Ok(())
}