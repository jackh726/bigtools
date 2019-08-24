use std::io::{self, Read, Seek, SeekFrom};
use std::io::{BufReader};
use std::fs::File;
use std::vec::Vec;

use byteordered::{ByteOrdered, Endianness};
use flate2::read::ZlibDecoder;

use crate::bigwig::{Value, ZoomHeader, CHROM_TREE_MAGIC, CIR_TREE_MAGIC, BIGWIG_MAGIC_LTH, BIGWIG_MAGIC_HTL};

/*
TODO
struct BigWigInterval<'a> {
    bigwig: &'a BigWigRead,
}

impl<'a> Iterator for BigWigInterval<'a> {
    type Item = io::Result<Value>;

    fn next(&mut self) -> Option<Self::Item> {
        None
    }
}
*/

#[derive(Debug)]
pub struct Block {
    pub offset: u64,
    pub size: u64,
}

#[derive(Clone, Debug)]
pub struct BBIHeader {
    pub endianness: Endianness,

    version: u16,
    zoom_levels: u16,
    chromosome_tree_offset: u64,
    full_data_offset: u64,
    full_index_offset: u64,
    field_count: u16,
    defined_field_count: u16,
    auto_sql_offset: u64,
    total_summary_offset: u64,
    uncompress_buf_size: u32,
}

#[derive(Clone, Debug)]
pub struct ChromInfo {
    pub name: String,
    id: u32,
    pub length: u32,
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
pub struct BigWigInfo {
    pub header: BBIHeader,
    zoom_headers: Vec<ZoomHeader>,
    chrom_info: Vec<ChromInfo>,
}

struct IntervalIter<'a> {
    bigwig: &'a mut BigWigRead,
    known_offset: u64,
    blocks: Vec<Block>,
    current_block: usize,
    vals: Option<Box<dyn Iterator<Item=Value> + Send + 'a>>,
}

impl<'a> Iterator for IntervalIter<'a> {
    type Item = io::Result<Value>;

    fn next(&mut self) -> Option<Self::Item> {
        let endianness = self.bigwig.info.header.endianness;
        let uncompress_buf_size: usize = self.bigwig.info.header.uncompress_buf_size as usize;

        let mut file = self.bigwig.reader.as_mut().unwrap();
        loop {
            match &mut self.vals {
                Some(vals) => {
                    match vals.next() {
                        Some(v) => {
                            return Some(Ok(v));
                        },
                        None => {
                            self.vals = None;
                            continue;
                        },
                    }
                },
                None => {
                    if self.current_block >= self.blocks.len() {
                        return None;
                    }
                    // TODO: Could minimize this by chunking block reads
                    let current_block = self.blocks.get(self.current_block).unwrap();
                    if self.known_offset != current_block.offset {
                        match file.seek(SeekFrom::Start(current_block.offset)) {
                            Ok(_) => {},
                            Err(e) => return Some(Err(e)),
                        }
                    }
                    let block_vals = match BigWigRead::get_block_values(&mut file, &current_block, endianness, uncompress_buf_size) {
                        Ok(vals) => vals,
                        Err(e) => return Some(Err(e)),
                    };
                    self.vals = Some(Box::new(block_vals));
                    self.known_offset = current_block.offset + current_block.size;
                    self.current_block += 1;
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

impl From<BigWigReadInfoError> for BigWigReadAttachError {
    fn from(error: BigWigReadInfoError) -> Self {
        return match error {
            BigWigReadInfoError::NotABigWig => BigWigReadAttachError::NotABigWig,
            BigWigReadInfoError::InvalidChroms => BigWigReadAttachError::InvalidChroms,
            BigWigReadInfoError::IoError(e) => BigWigReadAttachError::IoError(e),
        }
    }
}

pub enum BigWigReadInfoError {
    NotABigWig,
    InvalidChroms,
    IoError(io::Error),
}

impl From<io::Error> for BigWigReadInfoError {
    fn from(error: io::Error) -> Self {
        BigWigReadInfoError::IoError(error)
    }
}

pub struct BigWigRead {
    pub path: String,
    pub info: BigWigInfo,
    reader: Option<ByteOrdered<BufReader<File>, Endianness>>,
}

impl Clone for BigWigRead {
    fn clone(&self) -> Self {
        BigWigRead {
            path: self.path.clone(),
            info: self.info.clone(),
            reader: None,
        }
    }
}

impl BigWigRead {
    pub fn from_file_and_attach(path: String) -> Result<Self, BigWigReadAttachError> {
        let fp = File::open(path.clone())?;
        let file = BufReader::new(fp);
        let info = match BigWigRead::read_info(file) {
            Err(e) => {
                eprintln!("Error when opening: {}", path.clone());
                return Err(e.into());
            }
            Ok(info) => info,
        };

        Ok(BigWigRead {
            path,
            info,
            reader: None,
        })
    }

    pub fn get_chroms(&self) -> Vec<ChromAndSize> {
        self.info.chrom_info.iter().map(|c| ChromAndSize { name: c.name.clone(), length: c.length }).collect::<Vec<_>>()
    }

    pub fn ensure_reader(&mut self) -> io::Result<()> {
        if self.reader.is_none() {
            let endianness = self.info.header.endianness;
            let fp = File::open(self.path.clone())?;
            let file = ByteOrdered::runtime(BufReader::new(fp), endianness);
            self.reader.replace(file);
        }
        Ok(())
    }

    /// Manually close the open file descriptor (if it exists). If any operations are performed after this is called, the file descriptor will be reopened.
    pub fn close(&mut self) {
        self.reader.take();
    }

    fn read_info(file: BufReader<File>) -> Result<BigWigInfo, BigWigReadInfoError> {
        let mut file = ByteOrdered::runtime(file, Endianness::Little);

        let magic = file.read_u32()?;
        // println!("Magic {:x?}", magic);
        match magic {
            BIGWIG_MAGIC_HTL => file = file.into_opposite(),
            BIGWIG_MAGIC_LTH => {},
            _ => return Err(BigWigReadInfoError::NotABigWig),
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

        let zoom_headers = BigWigRead::read_zoom_headers(&mut file, &header)?;

        // TODO: could instead store this as an Option and only read when needed
        file.seek(SeekFrom::Start(header.chromosome_tree_offset))?;
        let magic = file.read_u32()?;
        let _block_size = file.read_u32()?;
        let key_size = file.read_u32()?;
        let val_size = file.read_u32()?;
        let item_count = file.read_u64()?;
        let _reserved = file.read_u64()?;
        if magic != CHROM_TREE_MAGIC {
            return Err(BigWigReadInfoError::InvalidChroms);
        }
        //println!("{:x?} {:?} {:?} {:?} {:?} {:?}", magic, _block_size, key_size, val_size, item_count, _reserved);
        assert_eq!(val_size, 8u32); 

        let mut chrom_info = Vec::with_capacity(item_count as usize);
        BigWigRead::read_chrom_tree_block(&mut file, &mut chrom_info, key_size)?;

        let info = BigWigInfo {
            header,
            zoom_headers,
            chrom_info,
        };

        Ok(info)
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
                BigWigRead::read_chrom_tree_block(f, chroms, key_size)?;
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
        BigWigRead::compare_position(chromq, chromq_start, chromb2, chromb2_end) <= 0 && BigWigRead::compare_position(chromq, chromq_end, chromb1, chromb1_start) >= 0
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
                if !BigWigRead::overlaps(chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base) {
                    continue;
                }
                //println!("Overlaps (leaf): {:?}:{:?}-{:?} with {:?}:{:?}-{:?}:{:?} {:?} {:?}", chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset, data_size);
                blocks.push(Block {
                    offset: data_offset,
                    size: data_size,
                })
            } else {
                let data_offset = file.read_u64()?;
                if !BigWigRead::overlaps(chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base) {
                    continue;
                }
                //println!("Overlaps (non-leaf): {:?}:{:?}-{:?} with {:?}:{:?}-{:?}:{:?} {:?}", chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset);
                childblocks.push(data_offset);
            }
        }
        for childblock in childblocks {
            //println!("Seeking to {:?}", childblock);
            file.seek(SeekFrom::Start(childblock))?;
            BigWigRead::search_overlapping_blocks(&mut file, chrom_ix, start, end, &mut blocks)?;
        }
        Ok(())
    }

    /// This assumes the file is at the cir tree start
    fn search_cir_tree(&mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<Vec<Block>> {
        let mut file = self.reader.as_mut().expect("reader must be seeked to cir tree start prior to calling search_cir_tree");

        let chrom_ix = {
            let chrom_info = &self.info.chrom_info;
            let chrom = chrom_info.iter().find(|&x| x.name == chrom_name);
            //println!("Chrom: {:?}", chrom);
            match chrom {
                Some(c) => c.id,
                None => return Err(std::io::Error::new(std::io::ErrorKind::Other, format!("{} not found.", chrom_name)))
            }
        };

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
        BigWigRead::search_overlapping_blocks(&mut file, chrom_ix, start, end, &mut blocks)?;
        //println!("overlapping_blocks: {:?}", blocks);
        Ok(blocks)
    }

    fn get_overlapping_blocks(&mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<Vec<Block>> {
        self.ensure_reader()?;
        let file = self.reader.as_mut().unwrap();

        let full_index_offset = self.info.header.full_index_offset;
        file.seek(SeekFrom::Start(full_index_offset))?;

        self.search_cir_tree(chrom_name, start, end)
    }

    /// This assumes that the file is currently at the block's start
    fn get_block_values(file: &mut ByteOrdered<BufReader<File>, Endianness>, block: &Block, endianness: Endianness, uncompress_buf_size: usize) -> io::Result<impl Iterator<Item=Value>> {
        let mut values: Vec<Value> = Vec::new();

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
        let _chrom_id = block_data_mut.read_u32()?;
        let chrom_start = block_data_mut.read_u32()?;
        let _chrom_end = block_data_mut.read_u32()?;
        let item_step = block_data_mut.read_u32()?;
        let item_span = block_data_mut.read_u32()?;
        let section_type = block_data_mut.read_u8()?;
        let _reserved = block_data_mut.read_u8()?;
        let item_count = block_data_mut.read_u16()?;

        let mut start = chrom_start;
        for _ in 0..item_count {
            match section_type {
                1 => {
                    // bedgraph
                    let chrom_start = block_data_mut.read_u32()?;
                    let chrom_end = block_data_mut.read_u32()?;
                    let value = block_data_mut.read_f32()?;
                    values.push(Value {
                        start: chrom_start,
                        end: chrom_end,
                        value,
                    });
                },
                2 => {
                    // variable step
                    let chrom_start = block_data_mut.read_u32()?;
                    let chrom_end = chrom_start + item_span;
                    let value = block_data_mut.read_f32()?;
                    values.push(Value {
                        start: chrom_start,
                        end: chrom_end,
                        value,
                    });
                },
                3 => {
                    // fixed step
                    let chrom_start = start;
                    start += item_step;
                    let chrom_end = chrom_start + item_span;
                    let value = block_data_mut.read_f32()?;
                    values.push(Value {
                        start: chrom_start,
                        end: chrom_end,
                        value,
                    });
                },
                _ => return Err(std::io::Error::new(std::io::ErrorKind::Other, format!("Unknown bigwig section type: {}", section_type)))
            }
        }

        Ok(values.into_iter())
    }

    // TODO: having problems with figuring out to use get_interval in bigwigmerge
    // This function only differs in three places:
    // 1) No 'a
    // 2) &'a mut self => mut self
    // 3) self.reader.as_mut().unwrap() => self.reader.unwrap()
    pub fn get_interval_move(mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<impl Iterator<Item=io::Result<Value>> + Send> {
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;

        let endianness = self.info.header.endianness;
        let uncompress_buf_size: usize = self.info.header.uncompress_buf_size as usize;

        self.ensure_reader()?;
        let mut file = self.reader.unwrap();

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
            let mut vals = || -> io::Result<Box<dyn Iterator<Item=io::Result<Value>> + Send>> {
                // TODO: Could minimize this by chunking block reads
                let vals = BigWigRead::get_block_values(&mut file, &block, endianness, uncompress_buf_size)?;
                if let Some(next_offset) = next_offset {
                    if next_offset != block.offset + block.size {
                        file.seek(SeekFrom::Start(next_offset))?;
                    }
                }
                Ok(Box::new(vals.map(|v| Ok(v))))
            };
            let v: Box<dyn Iterator<Item=io::Result<Value>> + Send> = vals().unwrap_or_else(|e| Box::new(std::iter::once(Err(e))));
            v
        }).filter_map(move |mut val| {
            if let Ok(ref mut v) = val {
                if v.end < start || v.start > end {
                    return None;
                }
                if v.start < start {
                    v.start = start;
                }
                if v.end > end {
                    v.end = end;
                }
            }
            return Some(val);
        });

        Ok(vals_iter)
    }

    pub fn get_interval<'a>(&'a mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<impl Iterator<Item=io::Result<Value>> + Send + 'a> {
        let blocks = self.get_overlapping_blocks(chrom_name, start, end)?;
        Ok(IntervalIter {
            bigwig: self,
            known_offset: 0,
            blocks: blocks,
            current_block: 0,
            vals: None,
        })
/*
        let endianness = self.info.header.endianness;
        let uncompress_buf_size: usize = self.info.header.uncompress_buf_size as usize;

        self.ensure_reader()?;
        let file = self.reader.as_mut().unwrap();

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
            let mut vals = || -> io::Result<Box<dyn Iterator<Item=io::Result<Value>> + Send + 'a>> {
                // TODO: Could minimize this by chunking block reads
                let vals = BigWigRead::get_block_values(file, &block, endianness, uncompress_buf_size)?;
                if let Some(next_offset) = next_offset {
                    if next_offset != block.offset + block.size {
                        file.seek(SeekFrom::Start(next_offset))?;
                    }
                }
                Ok(Box::new(vals.map(|v| Ok(v))))
            };
            let v: Box<dyn Iterator<Item=io::Result<Value>> + Send + 'a> = vals().unwrap_or_else(|e| Box::new(std::iter::once(Err(e))));
            v
        }).filter_map(move |mut val| {
            if let Ok(ref mut v) = val {
                if v.end <= start || v.start > end {
                    return None;
                }
                if v.start < start {
                    v.start = start;
                }
                if v.end > end {
                    v.end = end;
                }
            }
            return Some(val);
        });

        Ok(vals_iter)
*/
    }
}