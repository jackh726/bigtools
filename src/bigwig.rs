#![allow(non_snake_case)]
#![allow(dead_code)]

use std::io::{self, Read, Seek, SeekFrom, Write};
use std::io::{BufReader, BufWriter};
use std::fs::File;
use std::vec::Vec;

use futures::future::{Future, FutureExt, RemoteHandle};
use futures::channel::mpsc::{channel, Receiver};
use futures::executor::{block_on, LocalPool, ThreadPool};
use futures::sink::SinkExt;
use futures::stream::StreamExt;
use futures::task::SpawnExt;

use byteordered::{ByteOrdered, Endianness};

use byteorder::{NativeEndian, WriteBytesExt};

use flate2::Compression;
use flate2::write::ZlibEncoder;
use flate2::read::ZlibDecoder;

use serde::{Serialize, Deserialize};

use crate::chromvalues::ChromValues;
use crate::filebufferedchannel;
use crate::idmap::IdMap;
use crate::tell::Tell;
use crate::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};

const BIGWIG_MAGIC_LTH: u32 = 0x888F_FC26;
const BIGWIG_MAGIC_HTL: u32 = 0x26FC_8F88;
const BIGBED_MAGIC_LTH: u32 = 0x8789_F2EB;
const BIGBED_MAGIC_HTL: u32 = 0xEBF2_8987;

const CIR_TREE_MAGIC: u32 = 0x2468_ACE0;
const CHROM_TREE_MAGIC: u32 = 0x78CA_8C91;

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
    reserved: u64,
}

#[derive(Clone, Debug)]
struct ZoomHeader {
    reduction_level: u32,
    data_offset: u64,
    index_offset: u64,
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

#[derive(Debug)]
pub struct ChromSize {
    pub name: String,
    pub length: u32,
}

#[derive(Clone, Debug)]
pub struct BigWigInfo {
    pub header: BBIHeader,
    zoom_headers: Vec<ZoomHeader>,
    chrom_info: Vec<ChromInfo>,
}

#[derive(Debug)]
pub struct Block {
    pub offset: u64,
    pub size: u64,
}

#[derive(Debug, Clone)]
pub struct Value {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

#[derive(Debug, Clone)]
pub struct ValueWithChrom {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

#[derive(Debug)]
struct RTreeNode {
    start_chrom_idx: u32,
    start_base: u32,
    end_chrom_idx: u32,
    end_base: u32,
    children: RTreeChildren,
}

#[derive(Debug)]
enum RTreeChildren {
    DataSections(Vec<Section>),
    Nodes(Vec<RTreeNode>),
}

#[derive(Debug)]
struct SectionData {
    chrom: u32,
    start: u32,
    end: u32,
    data: Vec<u8>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Section {
    chrom: u32,
    start: u32,
    end: u32,
    offset: u64,
    size: u64,
}

#[derive(Debug)]
pub struct Summary {
    bases_covered: u64,
    min_val: f64,
    max_val: f64,
    sum: f64,
    sum_squares: f64,
}

type TempZoomInfo = (u32 /* resolution */, RemoteHandle<io::Result<usize>> /* Temp file that contains data */, TempFileBuffer, filebufferedchannel::Receiver<Section> /* sections */);
type ZoomInfo = (u32 /* resolution */, File /* Temp file that contains data */, Box<Iterator<Item=Section>> /* sections */);
const DEFAULT_ZOOM_SIZES: [u32; 11] = [10, 40, 160, 640, 2_560, 10_240, 40_960, 163_840, 655_360, 2_621_440, 10_485_760];

pub type ChromGroupRead = (
    Box<Future<Output=io::Result<Summary>> + Send + Unpin>,
    filebufferedchannel::Receiver<Section>,
    TempFileBuffer,
    Box<Future<Output=io::Result<usize>> + Send + Unpin>,
    Vec<TempZoomInfo>,
    (String, u32)
);

pub trait ChromGroupReadStreamingIterator {
    fn next(&mut self) -> io::Result<Option<ChromGroupRead>>;
}

#[derive(Debug)]
struct ZoomRecord {
    chrom: u32,
    start: u32,
    end: u32,
    valid_count: u32,
    min_value: f32,
    max_value: f32,
    sum: f32,
    sum_squares: f32,
}

struct BigWigInterval<'a> {
    bigwig: &'a BigWigRead,
}

impl<'a> Iterator for BigWigInterval<'a> {
    type Item = io::Result<Value>;

    fn next(&mut self) -> Option<Self::Item> {
        None
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
    pub fn from_file_and_attach(path: String) -> io::Result<Self> {
        let fp = File::open(path.clone())?;
        let file = BufReader::new(fp);
        let info = BigWigRead::read_info(file)?;

        Ok(BigWigRead {
            path,
            info,
            reader: None,
        })
    }

    pub fn get_chroms(&self) -> Vec<ChromAndSize> {
        self.info.chrom_info.iter().map(|c| ChromAndSize { name: c.name.clone(), length: c.length }).collect::<Vec<_>>()
    }

    pub fn ensure_reader<'a>(&'a mut self) -> io::Result<()> {
        if self.reader.is_none() {
            let endianness = self.info.header.endianness;
            let fp = File::open(self.path.clone())?;
            let file = ByteOrdered::runtime(BufReader::new(fp), endianness);
            self.reader.replace(file);
        }
        Ok(())
    }

    #[allow(clippy::all)]
    pub fn test_read_zoom(&mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<()> {
        let fp = File::open(self.path.clone())?;
        let file = BufReader::new(fp);

        if self.info.zoom_headers.is_empty() {
            println!("No zooms. Skipping test read.");
            return Ok(())
        }

        let uncompress_buf_size = self.info.header.uncompress_buf_size;
        let index_offset = self.info.zoom_headers[0].index_offset;
        let endianness = self.info.header.endianness;
        let mut file = ByteOrdered::runtime(file, endianness);
        file.seek(SeekFrom::Start(index_offset))?;

        let blocks = self.search_cir_tree(chrom_name, start, end)?;

        println!("Number of zoom blocks: {:?}", blocks.len());

        'blocks: for block in blocks {
            println!("Block: {:?}", block);
            file.seek(SeekFrom::Start(block.offset))?;

            let mut raw_data = vec![0u8; block.size as usize];
            file.read_exact(&mut raw_data)?;
            let data = if uncompress_buf_size > 0 {
                let mut uncompressed_block_data = vec![0u8; uncompress_buf_size as usize];
                let mut d = ZlibDecoder::new(&raw_data[..]);
                let _ = d.read(&mut uncompressed_block_data)?;
                uncompressed_block_data
            } else {
                raw_data
            };
            let itemcount = data.len() / (4 * 8);
            assert!(data.len() % (4 * 8) == 0);
            let mut data_mut = ByteOrdered::runtime(&data[..], endianness);
            for _ in 0..itemcount {
                let _chrom_id = data_mut.read_u32()?;
                let _chrom_start = data_mut.read_u32()?;
                let _chrom_end = data_mut.read_u32()?;
                let _valid_count = data_mut.read_u32()?;
                let _min_val = data_mut.read_f32()?;
                let _max_val = data_mut.read_f32()?;
                let _sum_data = data_mut.read_f32()?;
                let _sum_squares = data_mut.read_f32()?;
                println!("First zoom data: {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}", _chrom_id, _chrom_start, _chrom_end, _valid_count, _min_val, _max_val, _sum_data, _sum_squares);
                break 'blocks;
            }
        }

        Ok(())
    }

    fn read_info(file: BufReader<File>) -> io::Result<BigWigInfo> {
        let mut file = ByteOrdered::runtime(file, Endianness::Little);

        let magic = file.read_u32()?;
        //println!("Magic {:x?}: ", magic);
        match magic {
            BIGWIG_MAGIC_HTL => {
                file = file.into_opposite();
                true
            },
            BIGWIG_MAGIC_LTH => false,
            _ => return Err(std::io::Error::new(std::io::ErrorKind::Other, "File not a big wig"))
        };

        let version = file.read_u16()?;
        let zoom_levels = file.read_u16()?;
        let chromosome_tree_offset = file.read_u64()?;
        let full_data_offset = file.read_u64()?;
        let full_index_offset = file.read_u64()?;
        let field_count = file.read_u16()?;
        let defined_field_count = file.read_u16()?;
        let auto_sql_offset = file.read_u64()?;
        let total_summary_offset = file.read_u64()?;
        let uncompress_buf_size = file.read_u32()?;
        let reserved = file.read_u64()?;

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
            reserved,
        };

        println!("Header: {:?}", header);

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
            return Err(std::io::Error::new(std::io::ErrorKind::Other, "Invalid file format: CHROM_TREE_MAGIC does not match."))
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
                    Err(_) => return Err(std::io::Error::new(std::io::ErrorKind::Other, "Invalid file format: Invalid utf-8 string.")),
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

    fn search_cir_tree(&mut self, chrom_name: &str, start: u32, end: u32) -> io::Result<Vec<Block>> {
        self.ensure_reader()?;
        let mut file = self.reader.as_mut().unwrap();

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
}

#[derive(Clone)]
pub struct BigWigWriteOptions {
    pub compress: bool,
    pub items_per_slot: u32,
    pub block_size: u32,
}

pub struct BigWigWrite {
    pub path: String,
    pub options: BigWigWriteOptions,
}

impl BigWigWrite {
    pub fn create_file(path: String) -> io::Result<Self> {
        Ok(BigWigWrite {
            path,
            options: BigWigWriteOptions {
                compress: true,
                items_per_slot: 1024,
                block_size: 256,
            }
        })
    }

    const MAX_ZOOM_LEVELS: usize = 10;

    pub fn write_groups<V: 'static>(&self, chrom_sizes: std::collections::HashMap<String, u32>, vals: V) -> io::Result<()> where V : ChromGroupReadStreamingIterator + Send {
        let fp = File::create(self.path.clone())?;
        let mut file = BufWriter::new(fp);

        BigWigWrite::write_blank_headers(&mut file)?;

        let total_summary_offset = file.tell()?;
        file.write_all(&[0; 40])?;

        let full_data_offset = file.tell()?;

        {
            // Total items
            // Unless we know the vals ahead of time, we can't estimate total sections ahead of time.
            // Even then simply doing "(vals.len() as u32 + ITEMS_PER_SLOT - 1) / ITEMS_PER_SLOT"
            // underestimates because sections are split by chrom too, not just size.
            // Skip for now, and come back when we write real header + summary.
            file.write_u32::<NativeEndian>(0)?;
        }

        let pre_data = file.tell()?;
        let mut current_offset = pre_data;
        let (chrom_ids, summary, mut file, raw_sections_iter, zoom_infos) = block_on(self.write_vals(vals, file)?)?;
        let sections_iter = raw_sections_iter.map(|mut section| {
            // TODO: this assumes that all the data is contiguous
            // This will fail if we ever space the sections in any way
            section.offset = current_offset;
            current_offset += section.size;
            section
        });
        let (nodes, levels, total_sections) = BigWigWrite::get_rtreeindex(sections_iter, &self.options);
        let data_size = file.tell()? - pre_data;
        println!("Data size: {:?}", data_size);
        println!("Sections: {:?}", total_sections);
        println!("Summary: {:?}", summary);

        // Since the chrom tree is read before the index, we put this before the full data index
        // Therefore, there is a higher likelihood that the udc file will only need one read for chrom tree + full data index
        // Putting the chrom tree before the data also has a higher likelihood of being included with the beginning headers,
        //  but requires us to know all the data ahead of time (when writing)
        let chrom_index_start = file.tell()?;
        BigWigWrite::write_chrom_tree(&mut file, chrom_sizes, &chrom_ids.get_map())?;

        let index_start = file.tell()?;
        BigWigWrite::write_rtreeindex(&mut file, nodes, levels, total_sections, &self.options)?;

        let mut zoom_entries: Vec<ZoomHeader> = vec![];
        BigWigWrite::write_zooms(&mut file, zoom_infos, &mut zoom_entries, data_size, &self.options)?;

        //println!("Zoom entries: {:?}", zoom_entries);
        let num_zooms = zoom_entries.len() as u16;
        println!("Zooms: {:?}", num_zooms);

        // We *could* actually check the the real max size, but let's just assume at it's as large as the largest possible value
        // In most cases, I think this is the true max size (unless there is only one section and its less than ITEMS_PER_SLOT in size)
        let uncompress_buf_size = if self.options.compress {
            self.options.items_per_slot * (1 + 1 + 2 + 4 + 4 + 4 + 4 + 8 + 8)
        } else {
            0
        };

        file.seek(SeekFrom::Start(0))?;
        file.write_u32::<NativeEndian>(BIGWIG_MAGIC_LTH)?;
        file.write_u16::<NativeEndian>(4)?;
        file.write_u16::<NativeEndian>(num_zooms)?;
        file.write_u64::<NativeEndian>(chrom_index_start)?;
        file.write_u64::<NativeEndian>(full_data_offset)?;
        file.write_u64::<NativeEndian>(index_start)?;
        file.write_u16::<NativeEndian>(0)?;
        file.write_u16::<NativeEndian>(0)?;
        file.write_u64::<NativeEndian>(0)?;
        file.write_u64::<NativeEndian>(total_summary_offset)?;
        file.write_u32::<NativeEndian>(uncompress_buf_size)?;
        file.write_u64::<NativeEndian>(0)?;

        assert!(file.seek(SeekFrom::Current(0))? == 64);

        for zoom_entry in zoom_entries {
            file.write_u32::<NativeEndian>(zoom_entry.reduction_level)?;
            file.write_u32::<NativeEndian>(0)?;
            file.write_u64::<NativeEndian>(zoom_entry.data_offset)?;
            file.write_u64::<NativeEndian>(zoom_entry.index_offset)?;
        }

        file.seek(SeekFrom::Start(total_summary_offset))?;
        file.write_u64::<NativeEndian>(summary.bases_covered)?;
        file.write_f64::<NativeEndian>(summary.min_val)?;
        file.write_f64::<NativeEndian>(summary.max_val)?;
        file.write_f64::<NativeEndian>(summary.sum)?;
        file.write_f64::<NativeEndian>(summary.sum_squares)?;


        file.write_u32::<NativeEndian>(total_sections as u32)?;
        file.seek(SeekFrom::End(0))?;
        file.write_u32::<NativeEndian>(BIGWIG_MAGIC_LTH)?;

        Ok(())
    }

    fn write_blank_headers(file: &mut BufWriter<File>) -> io::Result<()> {
        file.seek(SeekFrom::Start(0))?;
        // Common header
        file.write_all(&[0; 64])?;
        // Zoom levels
        file.write_all(&[0; BigWigWrite::MAX_ZOOM_LEVELS * 24])?;

        Ok(())
    }

    fn write_chrom_tree(file: &mut BufWriter<File>, chrom_sizes: std::collections::HashMap<String, u32>, chrom_ids: &std::collections::HashMap<String, u32>) -> io::Result<()> {
        let mut chroms: Vec<&String> = chrom_ids.keys().collect();
        chroms.sort();
        //println!("Used chroms {:?}", chroms);

        file.write_u32::<NativeEndian>(CHROM_TREE_MAGIC)?;
        let item_count = chroms.len() as u64;
        // TODO: for now, always just use the length of chroms (if less than 256). This means we don't have to implement writing non-leaf nodes for now...
        // TODO: make this configurable
        let block_size = std::cmp::max(256, item_count) as u32;
        file.write_u32::<NativeEndian>(block_size)?;
        let max_bytes = chroms.iter().map(|a| a.len() as u32).fold(0, u32::max);
        file.write_u32::<NativeEndian>(max_bytes)?;
        file.write_u32::<NativeEndian>(8)?; // size of Id (u32) + Size (u32)
        file.write_u64::<NativeEndian>(item_count)?;
        file.write_u64::<NativeEndian>(0)?; // Reserved

        // Assuming this is all one block right now
        // TODO: add non-leaf nodes and split blocks
        file.write_u8(1)?;
        file.write_u8(0)?;
        file.write_u16::<NativeEndian>(item_count as u16)?;
        for chrom in chroms {
            let key_bytes = &mut vec![0u8; max_bytes as usize];
            let chrom_bytes = chrom.as_bytes();
            key_bytes[..chrom_bytes.len()].copy_from_slice(chrom_bytes);
            file.write_all(key_bytes)?;
            let id = *chrom_ids.get(chrom).expect("Internal error. (Chrom not found).");
            file.write_u32::<NativeEndian>(id)?;
            let length = chrom_sizes.get(&chrom[..]);
            match length {
                None => panic!("Expected length for chrom: {}", chrom),
                Some(l) => {
                    file.write_u32::<NativeEndian>(*l)?;
                }
            }
        }
        Ok(())
    }

    async fn write_section(compress: bool, items_in_section: Vec<Value>, chromId: u32) -> io::Result<SectionData> {
        let mut bytes: Vec<u8> = vec![];

        let start = items_in_section[0].start;
        let end = items_in_section[items_in_section.len() - 1].end;
        bytes.write_u32::<NativeEndian>(chromId)?;
        bytes.write_u32::<NativeEndian>(start)?;
        bytes.write_u32::<NativeEndian>(end)?;
        bytes.write_u32::<NativeEndian>(0)?;
        bytes.write_u32::<NativeEndian>(0)?;
        bytes.write_u8(1)?;
        bytes.write_u8(0)?;
        bytes.write_u16::<NativeEndian>(items_in_section.len() as u16)?;

        for item in items_in_section.iter() {
            bytes.write_u32::<NativeEndian>(item.start)?;
            bytes.write_u32::<NativeEndian>(item.end)?;
            bytes.write_f32::<NativeEndian>(item.value)?;   
        }

        let out_bytes = if compress {
            let mut e = ZlibEncoder::new(Vec::with_capacity(bytes.len()), Compression::default());
            e.write_all(&bytes)?;
            e.finish()?
        } else {
            bytes
        };

        Ok(SectionData {
            chrom: chromId,
            start,
            end,
            data: out_bytes,
        })
    }

    async fn write_zoom_section(compress: bool, items_in_section: Vec<ZoomRecord>) -> io::Result<SectionData> {
        let mut bytes: Vec<u8> = vec![];

        let start = items_in_section[0].start;
        let end = items_in_section[items_in_section.len() - 1].end;

        let chrom = items_in_section[0].chrom;
        for item in items_in_section.iter() {
            bytes.write_u32::<NativeEndian>(item.chrom)?;
            bytes.write_u32::<NativeEndian>(item.start)?;
            bytes.write_u32::<NativeEndian>(item.end)?;
            bytes.write_u32::<NativeEndian>(item.valid_count)?;
            bytes.write_f32::<NativeEndian>(item.min_value)?;
            bytes.write_f32::<NativeEndian>(item.max_value)?;
            bytes.write_f32::<NativeEndian>(item.sum)?;
            bytes.write_f32::<NativeEndian>(item.sum_squares)?; 
        }

        let out_bytes = if compress {
            let mut e = ZlibEncoder::new(Vec::with_capacity(bytes.len()), Compression::default());
            e.write_all(&bytes)?;
            e.finish()?
        } else {
            bytes
        };

        Ok(SectionData {
            chrom,
            start,
            end,
            data: out_bytes,
        })
    }

    pub fn read_group<I: 'static>(chrom: String, chromId: u32, mut group: I, mut pool: ThreadPool, options: BigWigWriteOptions)
        -> io::Result<ChromGroupRead>
        where I: ChromValues + Send {
        let cloned_chrom = chrom.clone();

        let num_zooms = DEFAULT_ZOOM_SIZES.len();

        let (mut ftx, frx) = channel::<_>(100);

        async fn create_do_write<W: Write>(mut data_file: W, mut section_sender: filebufferedchannel::Sender<Section>, mut frx: Receiver<impl Future<Output=io::Result<SectionData>>>) -> io::Result<usize> {
            let mut current_offset = 0;
            let mut total = 0;
            while let Some(section_raw) = frx.next().await {
                let section: SectionData = section_raw.await?;
                total += 1;
                let size = section.data.len() as u64;
                data_file.write_all(&section.data)?;
                section_sender.send(Section {
                    chrom: section.chrom,
                    start: section.start,
                    end: section.end,
                    offset: current_offset,
                    size: size,
                }).expect("Couldn't send section.");
                current_offset += size;
            }
            Ok(total)
        };

        let (sections_future, buf, section_receiver) = {
            let (buf, write) = TempFileBuffer::new()?;
            let file = BufWriter::new(write);

            let (section_sender, section_receiver) = filebufferedchannel::channel::<Section>(200);
            let sections_future = create_do_write(file, section_sender, frx);
            (sections_future, buf, section_receiver)
        };

        let process_zooms = move |zoom_channel: Receiver<_>, size: u32| -> io::Result<_> {
            let (buf, write) = TempFileBuffer::new()?;
            let file = BufWriter::new(write);

            let (section_sender, section_receiver) = filebufferedchannel::channel::<Section>(200);
            let file_future = create_do_write(file, section_sender, zoom_channel);

            Ok((size, file_future, buf, section_receiver))
        };

        let processed_zooms: Result<Vec<_>, _> = DEFAULT_ZOOM_SIZES.iter().map(|size| -> io::Result<_> {
            let (ftx, frx) = channel::<_>(100);
            let f = process_zooms(frx, *size)?;
            Ok((f, ftx))
        }).collect();
        let (zooms_futures, mut zooms_channels): (Vec<_>, Vec<_>) = processed_zooms?.into_iter().unzip();

        
        let read_file = async move || -> io::Result<Summary> {
            struct ZoomItem {
                live_info: Option<ZoomRecord>,
                records: Vec<ZoomRecord>,
            }
            struct BedGraphSection {
                items: Vec<Value>,
                zoom_items: Vec<ZoomItem>
            }

            let mut summary: Option<Summary> = None;

            let mut state_val = BedGraphSection {
                items: Vec::with_capacity(options.items_per_slot as usize),
                zoom_items: (0..num_zooms).map(|_| ZoomItem {
                    live_info: None,
                    records: Vec::with_capacity(options.items_per_slot as usize)
                }).collect(),
            };
            while let Some(current_val) = group.next()? {
                // TODO: test this correctly fails
                // TODO: change these to not panic
                assert!(current_val.start <= current_val.end);
                match group.peek() {
                    None => (),
                    Some(next_val) => {
                        assert!(
                            current_val.end <= next_val.start,
                            "Input bedGraph has overlapping values on chromosome {} at {}-{} and {}-{}",
                            chrom,
                            current_val.start,
                            current_val.end,
                            next_val.start,
                            next_val.end,
                        );
                    }
                }

                match &mut summary {
                    None => {
                        summary = Some(Summary {
                            bases_covered: u64::from(current_val.end - current_val.start),
                            min_val: f64::from(current_val.value),
                            max_val: f64::from(current_val.value),
                            sum: f64::from(current_val.end - current_val.start) * f64::from(current_val.value),
                            sum_squares: f64::from(current_val.end - current_val.start) * f64::from(current_val.value * current_val.value),
                        })
                    },
                    Some(summary) => {
                        summary.bases_covered += u64::from(current_val.end - current_val.start);
                        summary.min_val = summary.min_val.min(f64::from(current_val.value));
                        summary.max_val = summary.max_val.max(f64::from(current_val.value));
                        summary.sum += f64::from(current_val.end - current_val.start) * f64::from(current_val.value);
                        summary.sum_squares += f64::from(current_val.end - current_val.start) * f64::from(current_val.value * current_val.value);
                    }
                }

                for (i, mut zoom_item) in state_val.zoom_items.iter_mut().enumerate() {
                    let mut add_start = current_val.start;
                    loop {
                        if add_start >= current_val.end {
                            break
                        }
                        match &mut zoom_item.live_info {
                            None => {
                                zoom_item.live_info = Some(ZoomRecord {
                                    chrom: chromId,
                                    start: add_start,
                                    end: add_start,
                                    valid_count: 0,
                                    min_value: current_val.value,
                                    max_value: current_val.value,
                                    sum: 0.0,
                                    sum_squares: 0.0,
                                });
                            },
                            Some(zoom2) => {
                                let next_end = zoom2.start + DEFAULT_ZOOM_SIZES[i];
                                // End of bases that we could add
                                let add_end = std::cmp::min(next_end, current_val.end);
                                // If the last zoom ends before this value starts, we don't add anything
                                if add_end >= add_start {
                                    let added_bases = add_end - add_start;                                
                                    zoom2.end = add_end;
                                    zoom2.valid_count += added_bases;
                                    zoom2.min_value = zoom2.min_value.min(current_val.value);
                                    zoom2.max_value = zoom2.max_value.max(current_val.value);
                                    zoom2.sum += added_bases as f32 * current_val.value;
                                    zoom2.sum_squares += added_bases as f32 * current_val.value * current_val.value;
                                }
                                // If we made it to the end of the zoom (whether it was because the zoom ended before this value started,
                                // or we added to the end of the zoom), then write this zooms to the current section
                                if add_end == next_end {
                                    zoom_item.records.push(ZoomRecord {
                                        chrom: chromId,
                                        start: zoom2.start,
                                        end: zoom2.end,
                                        valid_count: zoom2.valid_count,
                                        min_value: zoom2.min_value,
                                        max_value: zoom2.max_value,
                                        sum: zoom2.sum,
                                        sum_squares: zoom2.sum_squares,
                                    });
                                    zoom_item.live_info = None;
                                }
                                // Set where we would start for next time
                                add_start = std::cmp::max(add_end, current_val.start);
                                // Write section if full
                                debug_assert!(zoom_item.records.len() <= options.items_per_slot as usize);
                                if zoom_item.records.len() == options.items_per_slot as usize {
                                    let items = std::mem::replace(&mut zoom_item.records, vec![]);
                                    let handle = pool.spawn_with_handle(BigWigWrite::write_zoom_section(options.compress, items)).expect("Couldn't spawn.");
                                    zooms_channels[i].send(handle.boxed()).await.expect("Couln't send");
                                }
                            }
                        }
                    }
                }
                state_val.items.push(current_val);
                if state_val.items.len() >= options.items_per_slot as usize {
                    let items = std::mem::replace(&mut state_val.items, vec![]);
                    let handle = pool.spawn_with_handle(BigWigWrite::write_section(options.compress, items, chromId)).expect("Couldn't spawn.");
                    ftx.send(handle.boxed()).await.expect("Couldn't send");
                }
            }

            if !state_val.items.is_empty() {
                let handle = pool.spawn_with_handle(BigWigWrite::write_section(options.compress, state_val.items, chromId)).expect("Couldn't spawn.");
                ftx.send(handle.boxed()).await.expect("Couldn't send");
            }

            for (i, mut zoom_item) in state_val.zoom_items.into_iter().enumerate() {
                if let Some(zoom2) = zoom_item.live_info {
                    assert!(chromId == zoom2.chrom);
                    zoom_item.records.push(ZoomRecord {
                        chrom: chromId,
                        start: zoom2.start,
                        end: zoom2.end,
                        valid_count: zoom2.valid_count,
                        min_value: zoom2.min_value,
                        max_value: zoom2.max_value,
                        sum: zoom2.sum,
                        sum_squares: zoom2.sum_squares,
                    });
                    zoom_item.live_info = None;
                }
                if !zoom_item.records.is_empty() {
                    let items = zoom_item.records;
                    let handle = pool.spawn_with_handle(BigWigWrite::write_zoom_section(options.compress, items)).expect("Couldn't spawn.");
                    zooms_channels[i].send(handle.boxed()).await.expect("Couln't send");
                }
            }

            let summary_complete = match summary {
                None => Summary {
                    bases_covered: 0,
                    min_val: 0.0,
                    max_val: 0.0,
                    sum: 0.0,
                    sum_squares: 0.0,
                },
                Some(summary) => summary,
            };
            Ok(summary_complete)
        };

        let (sections_remote, sections_handle) = sections_future.remote_handle();
        let (zoom_infos, zoom_remotes): (Vec<_>, Vec<_>) = zooms_futures.into_iter().map(|(size, file_future, buf, section_iter)| {
            let (remote, handle) = file_future.remote_handle();
            ((size, handle, buf, section_iter), remote)
        }).unzip();

        std::thread::spawn(move || {
            let mut pool = LocalPool::new();
            for zoom_remote in zoom_remotes {
                pool.spawner().spawn(zoom_remote).expect("Couldn't spawn future.");
            }
            pool.spawner().spawn(sections_remote).expect("Couldn't spawn future.");
            pool.run()
        });

        let (f_remote, f_handle) = read_file().remote_handle();
        std::thread::spawn(move || {
            block_on(f_remote);
        });
        Ok((Box::new(f_handle), section_receiver, buf, Box::new(sections_handle), zoom_infos, (cloned_chrom, chromId)))    
    }

    fn write_vals<V: 'static>(
        &self,
        mut vals_iter: V,
        file: BufWriter<File>
    ) -> io::Result<impl Future<Output=io::Result<(IdMap<String>, Summary, BufWriter<File>, Box<Iterator<Item=Section>>, Vec<ZoomInfo>)>>> where V : ChromGroupReadStreamingIterator + Send {
        let mut section_iter: Box<Iterator<Item=Section>> = Box::new(std::iter::empty());

        let mut zooms: Vec<(u32, Box<dyn Iterator<Item=Section>>, TempFileBuffer, TempFileBufferWriter)> = DEFAULT_ZOOM_SIZES.iter().map(|size| -> io::Result<_> {
            let section_iter: Box<Iterator<Item=Section>> = Box::new(std::iter::empty());

            let (buf, write) = TempFileBuffer::new()?;
            Ok((*size, section_iter, buf, write))
        }).collect::<Result<_, _>>()?;

        let mut raw_file = file.into_inner()?;

        let read_file = async move || -> io::Result<(IdMap<String>, Summary, BufWriter<File>, Box<Iterator<Item=Section>>, Vec<ZoomInfo>)> {
            let mut summary: Option<Summary> = None;

            let mut chrom_ids = IdMap::new();
            while let Some((summary_future, sections_receiver, mut sections_buf, sections_future, zoom_infos, (chrom, chrom_id))) = vals_iter.next()? {
                let real_id = chrom_ids.get_id(chrom);
                assert_eq!(real_id, chrom_id);
                sections_buf.switch(raw_file)?;

                let chrom_summary = summary_future.await?;
                let num_sections = sections_future.await?;
                dbg!(num_sections);
                section_iter = Box::new(section_iter.chain(sections_receiver.into_iter()));
                raw_file = sections_buf.await_file();

                for (i, (_size, future, buf, zoom_sections_receiver)) in zoom_infos.into_iter().enumerate() {
                    let zoom = &mut zooms[i];
                    assert_eq!(zoom.0, _size);
                    let num_sections = future.await?;
                    let mut old = std::mem::replace(&mut zoom.1, Box::new(std::iter::empty()));
                    old = Box::new(old.chain(zoom_sections_receiver.into_iter().take(num_sections)));
                    std::mem::replace(&mut zoom.1, old);
                    buf.expect_closed_write(&mut zoom.3)?;
                }

                match &mut summary {
                    None => {
                        summary = Some(chrom_summary)
                    },
                    Some(summary) => {
                        summary.bases_covered += chrom_summary.bases_covered;
                        summary.min_val = summary.min_val.min(chrom_summary.min_val);
                        summary.max_val = summary.max_val.max(chrom_summary.max_val);
                        summary.sum += chrom_summary.sum;
                        summary.sum_squares += chrom_summary.sum_squares;
                    }
                }
            }

            let summary_complete = summary.unwrap_or(Summary {
                bases_covered: 0,
                min_val: 0.0,
                max_val: 0.0,
                sum: 0.0,
                sum_squares: 0.0,
            });

            let zoom_infos: Vec<_> = zooms.into_iter().map(|zoom| {
                drop(zoom.3);
                (zoom.0, zoom.2.await_raw(), zoom.1)
            }).collect();
            Ok((chrom_ids, summary_complete, BufWriter::new(raw_file), section_iter, zoom_infos))
        };

        Ok(read_file())
    }

    fn write_zooms<'a>(mut file: &'a mut BufWriter<File>, zooms: Vec<ZoomInfo>, zoom_entries: &'a mut Vec<ZoomHeader>, data_size: u64, options: &BigWigWriteOptions) -> io::Result<()> {
        let mut zoom_count = 0;
        let mut last_zoom_section_count = u64::max_value();
        for zoom in zooms {
            let mut zoom_file = zoom.1;
            let zoom_size = zoom_file.seek(SeekFrom::End(0))?;
            if zoom_size > (data_size / 2) {
                //println!("Skipping zoom {:?} because it's size ({:?}) is greater than the data_size/2 ({:?})", zoom.0, zoom.3, data_size/2);
                continue;
            }
            let zoom_data_offset = file.tell()?;

            let mut current_offset = zoom_data_offset;
            let sections_iter = zoom.2.map(|mut section| {
                // TODO: assumes contiguous, see note for primary data
                section.offset = current_offset;
                current_offset += section.size;
                section
            });

            let (nodes, levels, total_sections) = BigWigWrite::get_rtreeindex(sections_iter, options);
            if last_zoom_section_count <= total_sections {
                continue;
            }
            last_zoom_section_count = total_sections;

            zoom_file.seek(SeekFrom::Start(0))?;
            let mut buf_reader = BufReader::new(zoom_file);
            std::io::copy(&mut buf_reader, &mut file)?;
            let zoom_index_offset = file.tell()?;
            //println!("Zoom {:?}, data: {:?}, offset {:?}", zoom.0, zoom_data_offset, zoom_index_offset);
            assert_eq!(zoom_index_offset - zoom_data_offset, zoom_size);
            BigWigWrite::write_rtreeindex(&mut file, nodes, levels, total_sections, options)?;

            zoom_entries.push(ZoomHeader {
                reduction_level: zoom.0,
                data_offset: zoom_data_offset,
                index_offset: zoom_index_offset,
            });

            zoom_count += 1;
            if zoom_count >= BigWigWrite::MAX_ZOOM_LEVELS {
                break;
            }
        }

        Ok(())
    }

    // TODO: it would be cool to output as an iterator so we don't have to store the index in memory
    fn get_rtreeindex<S>(sections_stream: S, options: &BigWigWriteOptions) -> (RTreeChildren, usize, u64) where S : Iterator<Item=Section> {
        use itertools::Itertools;

        let block_size = options.block_size as usize;
        let mut total_sections = 0;

        let sections: Vec<_> = sections_stream.collect();
        let chunks = sections.into_iter().inspect(|_| total_sections += 1).chunks(block_size);
        let mut current_nodes: Vec<RTreeChildren> = chunks.into_iter().map(|chunk| {
            let current_chunk: Vec<_> = chunk.collect();
            RTreeChildren::DataSections(current_chunk)
        }).collect();
        let mut levels = 0;
        let nodes: RTreeChildren = loop {
            if current_nodes.len() == 1 {
                break current_nodes.pop().unwrap();
            }
            levels += 1;
            let chunks = current_nodes.into_iter().chunks(block_size);
            current_nodes = chunks.into_iter().map(|chunk| {
                RTreeChildren::Nodes(chunk.map(|c| {
                    match &c {
                        RTreeChildren::DataSections(sections) => RTreeNode {
                            start_chrom_idx: sections.first().unwrap().chrom,
                            start_base: sections.first().unwrap().start,
                            end_chrom_idx: sections.last().unwrap().chrom,
                            end_base: sections.last().unwrap().end,
                            children: c,
                        },
                        RTreeChildren::Nodes(children) => RTreeNode {
                            start_chrom_idx: children.first().unwrap().start_chrom_idx,
                            start_base: children.first().unwrap().start_base,
                            end_chrom_idx: children.last().unwrap().end_chrom_idx,
                            end_base: children.last().unwrap().end_base,
                            children: c,
                        },
                    }
                }).collect())
            }).collect()
        };

        //println!("Total sections: {:?}", total_sections);
        //println!("Nodes ({:?}): {:?}", nodes.nodes.len(), nodes);
        //println!("Levels: {:?}", levels);
        (nodes, levels, total_sections)
    }


    fn write_rtreeindex(
        file: &mut BufWriter<File>,
        nodes: RTreeChildren,
        levels: usize,
        section_count: u64,
        options: &BigWigWriteOptions
    ) -> io::Result<()> {
        const NODEHEADER_SIZE: u64 = 1 + 1 + 2;
        const NON_LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8;
        const LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8 + 8;

        let mut index_offsets: Vec<u64> = vec![0u64; levels as usize];

        fn calculate_offsets(mut index_offsets: &mut Vec<u64>, nodes: &RTreeChildren, level: usize) {
            match nodes {
                RTreeChildren::DataSections(_) => (),
                RTreeChildren::Nodes(children) => {
                    index_offsets[level - 1] += NODEHEADER_SIZE;
                    for child in children {
                        index_offsets[level - 1] += NON_LEAFNODE_SIZE;
                        calculate_offsets(&mut index_offsets, &child.children, level - 1);
                    }
                }
            }
        }

        calculate_offsets(&mut index_offsets, &nodes, levels);
        //println!("index Offsets: {:?}", index_offsets);

        fn write_tree(
            mut file: &mut BufWriter<File>,
            nodes: &RTreeChildren,
            curr_level: usize,
            dest_level: usize,
            childnode_offset: u64,
            options: &BigWigWriteOptions
        ) -> io::Result<u64> {
            let NON_LEAFNODE_FULL_BLOCK_SIZE: u64 = NODEHEADER_SIZE + NON_LEAFNODE_SIZE * options.block_size as u64;
            let LEAFNODE_FULL_BLOCK_SIZE: u64 = NODEHEADER_SIZE + LEAFNODE_SIZE * options.block_size as u64;
            debug_assert!(curr_level >= dest_level);
            let mut total_size = 0;
            if curr_level != dest_level {
                let mut next_offset_offset = 0;
                match nodes {
                    RTreeChildren::DataSections(_) => panic!("Datasections found at level: {:?}", curr_level),
                    RTreeChildren::Nodes(children) => {
                        for child in children {
                            next_offset_offset += write_tree(&mut file, &child.children, curr_level - 1, dest_level, childnode_offset + next_offset_offset, options)?;
                        }
                    }
                }
                total_size += next_offset_offset;
                return Ok(total_size)
            }
            let isleaf = match nodes {
                RTreeChildren::DataSections(_) => 1,
                RTreeChildren::Nodes(_) => 0,
            };

            //println!("Level: {:?}", curr_level);
            //println!("Writing {}. Isleaf: {} At: {}", "node", isleaf, file.seek(SeekFrom::Current(0))?);
            file.write_u8(isleaf)?;
            file.write_u8(0)?;
            match &nodes {
                RTreeChildren::DataSections(sections) => {
                    file.write_u16::<NativeEndian>(sections.len() as u16)?;
                    total_size += 4;
                    for section in sections {
                        file.write_u32::<NativeEndian>(section.chrom)?;
                        file.write_u32::<NativeEndian>(section.start)?;
                        file.write_u32::<NativeEndian>(section.chrom)?;
                        file.write_u32::<NativeEndian>(section.end)?;
                        file.write_u64::<NativeEndian>(section.offset)?;
                        file.write_u64::<NativeEndian>(section.size)?;
                        total_size += 32;
                    }
                },
                RTreeChildren::Nodes(children) => {
                    file.write_u16::<NativeEndian>(children.len() as u16)?;
                    total_size += 4;
                    let full_size = if (curr_level - 1) > 0 {
                        NON_LEAFNODE_FULL_BLOCK_SIZE
                    } else {
                        LEAFNODE_FULL_BLOCK_SIZE
                    };
                    for (idx, child) in children.iter().enumerate() {
                        let child_offset: u64 = childnode_offset + idx as u64 * full_size;
                        file.write_u32::<NativeEndian>(child.start_chrom_idx)?;
                        file.write_u32::<NativeEndian>(child.start_base)?;
                        file.write_u32::<NativeEndian>(child.end_chrom_idx)?;
                        file.write_u32::<NativeEndian>(child.end_base)?;
                        file.write_u64::<NativeEndian>(child_offset)?;
                        total_size += 24;
                    }
                }
            }

            Ok(total_size)
        }


        let end_of_data = file.seek(SeekFrom::Current(0))?;
        {
            //println!("cirTree header (write):\n bs: {:?}\n ic: {:?}\n sci: {:?}\n sb: {:?}\n eci: {:?}\n eb: {:?}\n efo: {:?}\n ips: {:?}\n r: {:?}", BLOCK_SIZE, section_count, nodes.nodes[0].start_chrom_idx, nodes.nodes[0].start_base, nodes.nodes[nodes.nodes.len() - 1].end_chrom_idx, nodes.nodes[nodes.nodes.len() - 1].end_base, end_of_data, ITEMS_PER_SLOT, 0);
            file.write_u32::<NativeEndian>(CIR_TREE_MAGIC)?;
            file.write_u32::<NativeEndian>(options.block_size)?;
            file.write_u64::<NativeEndian>(section_count)?;
            match &nodes {
                RTreeChildren::DataSections(sections) => {
                    file.write_u32::<NativeEndian>(sections.first().unwrap().chrom)?;
                    file.write_u32::<NativeEndian>(sections.first().unwrap().start)?;
                    file.write_u32::<NativeEndian>(sections.last().unwrap().chrom)?;
                    file.write_u32::<NativeEndian>(sections.last().unwrap().end)?;
                },
                RTreeChildren::Nodes(children) => {
                    file.write_u32::<NativeEndian>(children.first().unwrap().start_chrom_idx)?;
                    file.write_u32::<NativeEndian>(children.first().unwrap().start_base)?;
                    file.write_u32::<NativeEndian>(children.last().unwrap().end_chrom_idx)?;
                    file.write_u32::<NativeEndian>(children.last().unwrap().end_base)?;
                }
            }
            file.write_u64::<NativeEndian>(end_of_data)?;
            file.write_u32::<NativeEndian>(options.items_per_slot)?;
            file.write_u32::<NativeEndian>(0)?;
        }

        let mut next_offset = file.seek(SeekFrom::Current(0))?;
        let mut total_size = 0;
        //println!("Levels: {:?}", levels);
        //println!("Start of index: {}", next_offset);
        for level in (0..=levels).rev() {
            if level > 0 {
                next_offset += index_offsets[level - 1];
            }
            total_size += write_tree(file, &nodes, levels, level, next_offset, options)?;
            //println!("End of index level {}: {}", level, file.seek(SeekFrom::Current(0))?);
        }
        //println!("Total index size: {:?}", total_size);

        Ok(())
    }

}
