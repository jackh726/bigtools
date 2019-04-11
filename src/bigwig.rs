#![allow(non_snake_case)]
#![allow(dead_code)]

use std::io::{Seek, SeekFrom};
use std::io::BufWriter;
use std::io::prelude::*;
use std::fs::File;
use std::vec::Vec;

use futures::future::{self, FutureExt};
use futures::channel::mpsc::{unbounded, UnboundedSender};
use futures::stream::{self, Stream, StreamExt};
use futures::task::SpawnExt;

use byteordered::{ByteOrdered, Endianness};

use byteorder::{NativeEndian, WriteBytesExt};

use flate2::Compression;
use flate2::write::ZlibEncoder;
use flate2::read::ZlibDecoder;

use crate::idmap::IdMap;

const BIGWIG_MAGIC_LTH: u32 = 0x888FFC26;
const BIGWIG_MAGIC_HTL: u32 = 0x26FC8F88;
const BIGBED_MAGIC_LTH: u32 = 0x8789F2EB;
const BIGBED_MAGIC_HTL: u32 = 0xEBF28987;

const CIR_TREE_MAGIC: u32 = 0x2468ACE0;
const CHROM_TREE_MAGIC: u32 = 0x78CA8C91;

#[derive(Debug)]
struct BBIHeader {
    endianness: Endianness,

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

#[derive(Debug)]
struct ZoomHeader {
    reduction_level: u32,
    data_offset: u64,
    index_offset: u64,
}

#[derive(Debug)]
struct ChromInfo {
    name: String,
    id: u32,
    length: u32,
}

#[derive(Debug)]
pub struct ChromSize {
    pub name: String,
    pub length: u32,
}

#[derive(Debug)]
pub struct BigWigInfo {
    header: Box<BBIHeader>,
    zoom_headers: Box<Vec<ZoomHeader>>,
    chrom_info: Box<Vec<ChromInfo>>,
}

#[derive(Debug)]
struct Block {
    offset: u64,
    size: u64,
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
struct RTreeNodeList<RTreeNode> {
    nodes: Vec<RTreeNode>
}

#[derive(Debug)]
struct RTreeNode {
    start_chrom_idx: u32,
    start_base: u32,
    end_chrom_idx: u32,
    end_base: u32,
    kind: RTreeNodeType,
}

#[derive(Debug)]
enum RTreeNodeType {
    Leaf {
        offset: u64,
        size: u64,
    },
    NonLeaf {
        children: RTreeNodeList<RTreeNode>,
    },
}

#[derive(Debug)]
struct SectionData {
    chrom: u32,
    start: u32,
    end: u32,
    data: Vec<u8>,
}

#[derive(Debug)]
struct Section {
    offset: u64,
    size: u64,
    chrom: u32,
    start: u32,
    end: u32,
}

#[derive(Debug)]
struct Summary {
    bases_covered: u64,
    min_val: f64,
    max_val: f64,
    sum: f64,
    sum_squares: f64,
}

type TempZoomInfo = (u32 /* resolution */, File /* Temp file that contains data */, Vec<Section> /* sections */);

struct BedGraphSectionItem {
    start: u32,
    end: u32,
    val: f32,
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

pub struct BigWig {
    pub path: String,
    fp: File,
    info: Option<Box<BigWigInfo>>,
}

impl BigWig {
    pub fn test_read_zoom(&mut self, chrom_name: &str, start: u32, end: u32) -> std::io::Result<()> {
        if self.info.as_ref().unwrap().zoom_headers.len() == 0 {
            println!("No zooms. Skipping test read.");
            return Ok(())
        }
        let chrom_ix = {
            let info = self.info.as_ref().unwrap();
            let chrom_info = &info.chrom_info;
            let chrom = chrom_info.iter().find(|&x| x.name == chrom_name);
            println!("Chrom: {:?}", chrom);
            match chrom {
                Some(c) => c.id,
                None => return Err(std::io::Error::new(std::io::ErrorKind::Other, format!("{} not found.", chrom_name)))
            }
        };
        let uncompress_buf_size = self.info.as_ref().unwrap().header.uncompress_buf_size;
        let index_offset = self.info.as_ref().unwrap().zoom_headers[0].index_offset;
        let endianness = self.info.as_ref().unwrap().header.endianness;
        let mut file = ByteOrdered::runtime(self.get_buf_reader(), endianness);
        file.seek(SeekFrom::Start(index_offset))?;
        let blocks = BigWig::get_overlapping_blocks(&mut file, chrom_ix, start, end)?;

        for block in blocks {
            file.seek(SeekFrom::Start(block.offset))?;

            let mut raw_data = vec![0u8; block.size as usize];
            file.read_exact(&mut raw_data)?;
            let data = if uncompress_buf_size > 0 {
                let mut uncompressed_block_data = vec![0u8; uncompress_buf_size as usize];
                let mut d = ZlibDecoder::new(&raw_data[..]);
                d.read(&mut uncompressed_block_data)?;
                uncompressed_block_data
            } else {
                raw_data
            };
            let itemcount = data.len() / (4 * 8);
            assert!(data.len() % (4 * 8) == 0);
            let mut data_mut = ByteOrdered::runtime(&data[..], endianness);
            for _ in 0..itemcount {
                data_mut.read_u32()?;
                data_mut.read_u32()?;
                data_mut.read_u32()?;
                data_mut.read_u32()?;
                data_mut.read_f32()?;
                data_mut.read_f32()?;
                data_mut.read_f32()?;
                data_mut.read_f32()?;
            }
        }

        Ok(())
    }
    pub fn from_file(path: String) -> std::io::Result<Self> {
        let fp = File::open(path.clone())?;
        println!("{:?}", fp);
        return Ok(
            BigWig {
                path,
                fp: fp,
                info: Option::None,
            }
        )
    }

    pub fn create_file(path: String) -> std::io::Result<Self> {
        let fp = File::create(path.clone())?;
        println!("{:?}", fp);
        return Ok(
            BigWig {
                path,
                fp: fp,
                info: Option::None,
            }
        )
    }

    pub fn read_info(&mut self) -> std::io::Result<()> {
        if let Some(_) = &self.info {
            return Ok(());
        }

        let mut file = ByteOrdered::runtime(self.get_buf_reader(), Endianness::Little);

        let magic = file.read_u32()?;
        println!("Magic {:x?}: ", magic);
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

        let zoom_headers = BigWig::read_zoom_headers(&mut file, &header)?;

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
        BigWig::read_chrom_tree_block(&mut file, &mut chrom_info, key_size)?;

        let info = BigWigInfo {
            header: Box::new(header),
            zoom_headers: Box::new(zoom_headers),
            chrom_info: Box::new(chrom_info),
        };

        self.info = Some(Box::new(info));

        println!("Info read successfully.");
        Ok(())
    }

    fn ensure_info(&self) -> Result<(), &'static str> {
        if let Some(_) = &self.info {
            return Ok(());
        }
        Err("Must first call read_info")
    }

    fn read_zoom_headers(file: &mut ByteOrdered<std::io::BufReader<&File>, Endianness>, header: &BBIHeader) -> std::io::Result<Vec<ZoomHeader>> {
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

    fn read_chrom_tree_block(f: &mut ByteOrdered<std::io::BufReader<&File>, Endianness>, chroms: &mut Vec<ChromInfo>, key_size: u32) -> std::io::Result<()> {
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
                BigWig::read_chrom_tree_block(f, chroms, key_size)?;
                f.seek(SeekFrom::Start(current_position))?;
            }
        }
        return Ok(())
    }

    #[inline]
    fn compare_position(chrom1: u32, chrom1_base: u32, chrom2: u32, chrom2_base: u32) -> i8 {
        if chrom1 < chrom2 {
            return -1;
        } else if chrom1 > chrom2 {
            return 1;
        } else {
            if chrom1_base < chrom2_base {
                return -1;
            } else if chrom1_base > chrom2_base {
                return 1;
            } else {
                return 0;
            }
        }
    }

    fn overlaps(chromq: u32, chromq_start: u32, chromq_end: u32, chromb1: u32, chromb1_start: u32, chromb2: u32, chromb2_end: u32) -> bool {
        return BigWig::compare_position(chromq, chromq_start, chromb2, chromb2_end) <= 0 && BigWig::compare_position(chromq, chromq_end, chromb1, chromb1_start) >= 0;
    }

    fn search_overlapping_blocks(mut file: &mut ByteOrdered<std::io::BufReader<&File>, Endianness>, chrom_ix: u32, start: u32, end: u32, mut blocks: &mut Vec<Block>) -> std::io::Result<()> {
        //println!("Searching for overlapping blocks at {:?}. Searching {:?}:{:?}-{:?}", self.current_file_offset()?, chrom_ix, start, end);

        let isleaf: u8 = file.read_u8()?;
        assert!(isleaf == 1 || isleaf == 0, "Unexpected isleaf: {}", isleaf);
        let _reserved = file.read_u8()?;
        let count: u16 = file.read_u16()?;
        println!("Index: {:?} {:?} {:?}", isleaf, _reserved, count);

        let mut childblocks: Vec<u64> = vec![];
        for _ in 0..count {
            let start_chrom_ix = file.read_u32()?;
            let start_base = file.read_u32()?;
            let end_chrom_ix = file.read_u32()?;
            let end_base = file.read_u32()?;
            if isleaf == 1 {
                let data_offset = file.read_u64()?;
                let data_size = file.read_u64()?;
                if !BigWig::overlaps(chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base) {
                    continue;
                }
                println!("Overlaps (leaf): {:?}:{:?}-{:?} with {:?}:{:?}-{:?}:{:?} {:?} {:?}", chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset, data_size);
                blocks.push(Block {
                    offset: data_offset,
                    size: data_size,
                })
            } else {
                let data_offset = file.read_u64()?;
                if !BigWig::overlaps(chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base) {
                    continue;
                }
                println!("Overlaps (non-leaf): {:?}:{:?}-{:?} with {:?}:{:?}-{:?}:{:?} {:?}", chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset);
                childblocks.push(data_offset);
            }
        }
        for childblock in childblocks {
            println!("Seeking to {:?}", childblock);
            file.seek(SeekFrom::Start(childblock))?;
            BigWig::search_overlapping_blocks(&mut file, chrom_ix, start, end, &mut blocks)?;
        }
        return Ok(());
    }

    fn get_overlapping_blocks(file: &mut ByteOrdered<std::io::BufReader<&File>, Endianness>, chrom_ix: u32, start: u32, end: u32) -> std::io::Result<Vec<Block>> {
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

        println!("cirTree header:\n bs: {:?}\n ic: {:?}\n sci: {:?}\n sb: {:?}\n eci: {:?}\n eb: {:?}\n efo: {:?}\n ips: {:?}\n r: {:?}", _blocksize, _item_count, _start_chrom_idx, _start_base, _end_chrom_idx, _end_base, _end_file_offset, _item_per_slot, _reserved);
        let mut blocks: Vec<Block> = vec![];
        BigWig::search_overlapping_blocks(file, chrom_ix, start, end, &mut blocks)?;
        println!("overlapping_blocks: {:?}", blocks);
        Ok(blocks)
    }

    pub fn get_interval(&mut self, chrom_name: &str, start: u32, end: u32) -> std::io::Result<Vec<Value>> {
        self.ensure_info().or(Err(std::io::Error::new(std::io::ErrorKind::Other, "Must first call read_info")))?;

        let chrom_ix = {
            let info = self.info.as_ref().unwrap();
            let chrom_info = &info.chrom_info;
            let chrom = chrom_info.iter().find(|&x| x.name == chrom_name);
            println!("Chrom: {:?}", chrom);
            match chrom {
                Some(c) => c.id,
                None => return Err(std::io::Error::new(std::io::ErrorKind::Other, format!("{} not found.", chrom_name)))
            }
        };
        let uncompress_buf_size: usize = {
            let info = self.info.as_ref().unwrap();
            info.header.uncompress_buf_size as usize
        };
        let full_index_offset = self.info.as_ref().unwrap().header.full_index_offset;
        let endianness = self.info.as_ref().unwrap().header.endianness;
        let mut file = ByteOrdered::runtime(self.get_buf_reader(), endianness);
        file.seek(SeekFrom::Start(full_index_offset))?;
        let blocks = BigWig::get_overlapping_blocks(&mut file, chrom_ix, start, end)?;

        let mut values: Vec<Value> = Vec::new();

        for block in blocks {
            // TODO: Could minimize this by chunking block reads
            file.seek(SeekFrom::Start(block.offset))?;
            let mut raw_data = vec![0u8; block.size as usize];
            file.read_exact(&mut raw_data)?;
            let block_data: Vec<u8> = if uncompress_buf_size > 0 {
                let mut uncompressed_block_data = vec![0u8; uncompress_buf_size];
                let mut d = ZlibDecoder::new(&raw_data[..]);
                d.read(&mut uncompressed_block_data)?;
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
                            start: std::cmp::max(chrom_start, start),
                            end: std::cmp::min(chrom_end, end),
                            value,
                        });
                    },
                    2 => {
                        // variable step
                        let chrom_start = block_data_mut.read_u32()?;
                        let chrom_end = chrom_start + item_span;
                        let value = block_data_mut.read_f32()?;
                        values.push(Value {
                            start: std::cmp::max(chrom_start, start),
                            end: std::cmp::min(chrom_end, end),
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
                            start: std::cmp::max(chrom_start, start),
                            end: std::cmp::min(chrom_end, end),
                            value,
                        });
                    },
                    _ => return Err(std::io::Error::new(std::io::ErrorKind::Other, format!("Unknown bigwig section type: {}", section_type)))
                }
            }
        }
        Ok(values)
    }

    const MAX_ZOOM_LEVELS: usize = 10;

    pub fn write<V>(&mut self, chrom_sizes: std::collections::HashMap<String, u32>, mut vals: V) -> std::io::Result<()> where V : std::iter::Iterator<Item=ValueWithChrom> + std::marker::Send {
        self.write_blank_headers()?;

        let total_summary_offset = self.current_file_offset()?;

        self.write_blank_summary()?;

        let chrom_ids = IdMap::new();

        let full_data_offset = self.current_file_offset()?;

        {
            // Total items
            // Unless we know the vals ahead of time, we can't estimate total sections ahead of time.
            // Even then simply doing "(vals.len() as u32 + ITEMS_PER_SLOT - 1) / ITEMS_PER_SLOT"
            // underestimates because sections are split by chrom too, not just size.
            // Skip for now, and come back when we write real header + summary.
            (&mut self.fp).write_u32::<NativeEndian>(0)?;
        }

        let pre_data = self.current_file_offset()?;
        // Ideally, should instead be passing a mutable borrow to self.fp
        // However, write_vals needs ownership because of threading
        // Alternatively, could not have a self.fp field and instead pass file around explicitly
        let fp = std::mem::replace(&mut self.fp, tempfile::tempfile()?);
        let (chrom_summary_future, sections_stream, zooms_future, file_future) = self.write_vals(&mut vals, chrom_ids, fp)?;
        let rti_future = BigWig::get_rtreeindex(sections_stream);
        let fs = async {
            futures::join!(chrom_summary_future, rti_future)
        };
        let (cs, rti) = futures::executor::block_on(fs);
        let (chrom_ids, summary) = cs;
        let (nodes, levels, total_sections) = rti;
        let zooms = futures::executor::block_on(zooms_future);
        let file = futures::executor::block_on(file_future);
        std::mem::replace(&mut self.fp, file);
        let data_size = self.current_file_offset()? - pre_data;
        println!("Data size: {:?}", data_size);
        println!("Sections: {:?}", total_sections);
        println!("Summary: {:?}", summary);
        println!("Zooms: {:?}", zooms.len());

        // Since the chrom tree is read before the index, we put this before the full data index
        // Therefore, there is a higher likelihood that the udc file will only need one read for chrom tree + full data index
        // Putting the chrom tree before the data also has a higher likelihood of being included with the beginning headers,
        //  but requires us to know all the data ahead of time (when writing)
        let chrom_index_start = self.current_file_offset()?;
        self.write_chrom_tree(chrom_sizes, &chrom_ids.get_map())?;

        let index_start = self.current_file_offset()?;
        self.write_rtreeindex(nodes, levels, total_sections)?;

        const DO_COMPRESS: bool = true; // TODO: param
        const ITEMS_PER_SLOT: u32 = 1024; // TODO: param

        let mut zoom_entries: Vec<ZoomHeader> = vec![];
        self.write_zooms(zooms, &mut zoom_entries)?;

        //println!("Zoom entries: {:?}", zoom_entries);
        let num_zooms = zoom_entries.len() as u16;

        // We *could* actually check the the real max size, but let's just assume at it's as large as the largest possible value
        // In most cases, I think this is the true max size (unless there is only one section and its less than ITEMS_PER_SLOT in size)
        let uncompress_buf_size = if DO_COMPRESS {
            ITEMS_PER_SLOT * (1 + 1 + 2 + 4 + 4 + 4 + 4 + 8 + 8)
        } else {
            0
        };

        let mut file = self.get_buf_writer();
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

    fn current_file_offset(&mut self) -> std::io::Result<u64> {
        let file = &mut self.fp;
        file.seek(SeekFrom::Current(0))
    }

    fn get_buf_reader(&mut self) -> std::io::BufReader<&File> {
        let file = &mut self.fp;
        std::io::BufReader::new(file)
    }

    fn get_buf_writer(&mut self) -> std::io::BufWriter<&File> {
        let file = &mut self.fp;
        //std::io::BufWriter::new(file)
        std::io::BufWriter::with_capacity(1024 * 1024 * 5, file)
    }

    fn write_blank_headers(&mut self) -> std::io::Result<()> {
        let mut file = self.get_buf_writer();
        file.seek(SeekFrom::Start(0))?;
        // Common header
        file.write_all(&vec![0; 64])?;
        // Zoom levels
        file.write_all(&vec![0; BigWig::MAX_ZOOM_LEVELS * 24])?;

        Ok(())
    }

    fn write_blank_summary(&mut self) -> std::io::Result<()> {
        let mut file = self.get_buf_writer();
        file.write_all(&vec![0; 40])?;

        Ok(())
    }

    fn write_chrom_tree(&mut self, chrom_sizes: std::collections::HashMap<String, u32>, chrom_ids: &std::collections::HashMap<String, u32>) -> std::io::Result<()> {
        let mut chroms: Vec<&String> = chrom_ids.keys().collect();
        chroms.sort();
        println!("Used chroms {:?}", chroms);

        let mut file = self.get_buf_writer();
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
            let id = *chrom_ids.get(chrom).unwrap();
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

    fn write_section(items_in_section: Vec<BedGraphSectionItem>, chromId: u32) -> std::io::Result<SectionData> {
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
            bytes.write_f32::<NativeEndian>(item.val)?;   
        }

        let COMPRESS = true;

        let out_bytes = if COMPRESS {
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

    fn write_zoom_section(items_in_section: Vec<ZoomRecord>, chromId: u32) -> std::io::Result<SectionData> {
        let mut bytes: Vec<u8> = vec![];

        let start = items_in_section[0].start;
        let end = items_in_section[items_in_section.len() - 1].end;

        for item in items_in_section.iter() {
            bytes.write_u32::<NativeEndian>(chromId)?;
            bytes.write_u32::<NativeEndian>(item.start)?;
            bytes.write_u32::<NativeEndian>(item.end)?;
            bytes.write_u32::<NativeEndian>(item.valid_count)?;
            bytes.write_f32::<NativeEndian>(item.min_value)?;
            bytes.write_f32::<NativeEndian>(item.max_value)?;
            bytes.write_f32::<NativeEndian>(item.sum)?;
            bytes.write_f32::<NativeEndian>(item.sum_squares)?; 
        }

        let COMPRESS = true;

        let out_bytes = if COMPRESS {
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

    fn write_vals<'a, V>(
        &mut self,
        vals_iter: V,
        mut chrom_ids: IdMap<String>,
        mut file: File
    ) -> std::io::Result<(
        impl futures::Future<Output=(IdMap<String>, Summary)>,
        impl futures::stream::Stream<Item=Section>,
        impl futures::Future<Output=Vec<TempZoomInfo>>,
        impl futures::Future<Output=File>
        )> where V : std::iter::Iterator<Item=ValueWithChrom> + std::marker::Send {
        const ITEMS_PER_SLOT: u16 = 1024;

        #[derive(Debug)]
        struct LiveZoomInfo {
            chrom: String,
            start: u32,
            end: u32,
            valid_count: u32,
            min_value: f32,
            max_value: f32,
            sum: f32,
            sum_squares: f32,
        }

        type ZoomInfo<'a> = (u32, BufWriter<File>, Option<LiveZoomInfo>);

        let zoom_sizes: Vec<u32> = vec![10, 40, 160, 640, 2560, 10240, 40960, 163840, 655360, 2621440];
        let mut zooms: Vec<ZoomInfo> = vec![];
        for size in zoom_sizes.iter() {
            let temp = tempfile::tempfile()?;
            zooms.push((*size, std::io::BufWriter::new(temp), None));
        }
        let num_zooms = zooms.len();
        let mut zoom_sections_out: Vec<Vec<Section>> = zoom_sizes.iter().map(|_| vec![]).collect();

        let (mut ftx, mut frx) = unbounded::<_>();
        let (ftxsections, frxsections) = unbounded::<_>();
        let (ftxzooms, mut frxzooms) = unbounded::<(u32, _)>();

        let do_write = async move || -> std::io::Result<File> {
            let mut current_offset = file.seek(SeekFrom::End(0))?;
            while let Some(section_raw) = await!(frx.next()) {
                let section: SectionData = await!(section_raw)?;
                let size = section.data.len() as u64;
                let section_offset = current_offset;
                current_offset += size;
                file.write_all(&section.data)?;
                ftxsections.unbounded_send(Section {
                    chrom: section.chrom,
                    start: section.start,
                    end: section.end,
                    offset: section_offset,
                    size,
                }).unwrap();
            }
            Ok(file)
        };


        let sections_future = do_write();
        let read_file = async move || -> std::io::Result<(IdMap<String>, Summary)> {
            let mut summary: Option<Summary> = None;
            
            struct ZoomItem {
                live_info: Option<LiveZoomInfo>,
                records: Vec<ZoomRecord>,
            }
            struct BedGraphSection {
                chrom: String,
                items: Vec<BedGraphSectionItem>,
                zoom_items: Vec<ZoomItem>
            }

            let pool = futures::executor::ThreadPoolBuilder::new().pool_size(4).create().expect("Unable to create thread pool.");
            let mut pool1 = pool.clone();
            let mut pool2 = pool.clone();

            let mut state: Option<BedGraphSection> = None;

            let mut write_section_items = |items: Vec<BedGraphSectionItem>, chromId: u32, ftx: &mut UnboundedSender<_>| {
                let (remote, handle) = future::lazy(move |_| {
                    BigWig::write_section(items, chromId)
                }).remote_handle();
                pool1.spawn(remote).expect("Couldn't spawn.");
                ftx.unbounded_send(handle.boxed()).expect("Couldn't send");
            };

            let mut write_zoom_items = |i: usize, items: Vec<ZoomRecord>| {
                let chromId = items[0].chrom;
                let (remote, handle) = future::lazy(move |_| {
                    BigWig::write_zoom_section(items, chromId)
                }).remote_handle();
                pool2.spawn(remote).expect("Couldn't spawn.");
                ftxzooms.unbounded_send((i as u32, handle.boxed())).expect("Couldn't send");
            };
            for current_val in vals_iter {
                // TODO: test this correctly fails
                if let Some(state) = &state {
                    assert!(state.chrom <= current_val.chrom, "Input bedGraph not sorted by chromosome. Sort with `sort -k1,1 -k2,2n`.");
                    if state.chrom == current_val.chrom {
                        let lastItem = state.items.last().expect("Invalid state");
                        assert!(
                            lastItem.end <= current_val.start,
                            "Input bedGraph has overlapping values on chromosome {} at {}-{} and {}-{}",
                            current_val.chrom,
                            lastItem.start,
                            lastItem.end,
                            current_val.start,
                            current_val.end,
                        );
                    }
                }
                assert!(current_val.start <= current_val.end);

                if let Some(state_val) = &mut state {
                    let diffchrom = current_val.chrom != state_val.chrom;
                    if diffchrom {
                        let mut state_val = state.unwrap();
                        let chromId = chrom_ids.get_id(state_val.chrom);
                        let items = std::mem::replace(&mut state_val.items, vec![]);
                        write_section_items(items, chromId, &mut ftx);
                        for (i, mut zoom_item) in state_val.zoom_items.into_iter().enumerate() {
                            if let Some(zoom2) = zoom_item.live_info {
                                let chromId = chrom_ids.get_id(zoom2.chrom.to_string());
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
                            write_zoom_items(i, zoom_item.records);
                        }
                        state = None;
                    }
                }
                if let None = state {
                    state = Some(BedGraphSection {
                        chrom: current_val.chrom.clone(),
                        items: Vec::with_capacity(ITEMS_PER_SLOT as usize),
                        zoom_items: (0..num_zooms).map(|_| ZoomItem {
                            live_info: None,
                            records: Vec::with_capacity(ITEMS_PER_SLOT as usize)
                        }).collect(),
                    });
                }

                let state_val = state.as_mut().unwrap();
                if state_val.items.len() >= ITEMS_PER_SLOT as usize {
                    let chromId = chrom_ids.get_id(state_val.chrom.clone());
                    let items = std::mem::replace(&mut state_val.items, vec![]);
                    write_section_items(items, chromId, &mut ftx);
                }
                for (i, mut zoom_item) in state_val.zoom_items.iter_mut().enumerate() {
                    let mut add_start = current_val.start;
                    loop {
                        if add_start >= current_val.end {
                            break
                        }
                        match &mut zoom_item.live_info {
                            None => {
                                zoom_item.live_info = Some(LiveZoomInfo {
                                    chrom: current_val.chrom.clone(),
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
                                let next_end = zoom2.start + zoom_sizes[i];
                                // End of bases that we could add
                                let add_end = std::cmp::min(next_end, current_val.end);
                                // If the last zoom ends before this value starts, we don't add anything
                                if add_end >= add_start {
                                    let added_bases = add_end - add_start;                                
                                    zoom2.end = current_val.end;
                                    zoom2.valid_count += added_bases;
                                    zoom2.min_value = zoom2.min_value.min(current_val.value);
                                    zoom2.max_value = zoom2.max_value.max(current_val.value);
                                    zoom2.sum += added_bases as f32 * current_val.value;
                                    zoom2.sum_squares += added_bases as f32 * current_val.value * current_val.value;
                                }
                                // If we made it to the end of the zoom (whether it was because the zoom ended before this value started,
                                // or we added to the end of the zoom), then write this zooms to the current section
                                if add_end == next_end {
                                    let chromId = chrom_ids.get_id(zoom2.chrom.to_string());
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
                                assert!(zoom_item.records.len() <= ITEMS_PER_SLOT as usize);
                                if zoom_item.records.len() == ITEMS_PER_SLOT as usize {
                                    let items = std::mem::replace(&mut zoom_item.records, vec![]);
                                    write_zoom_items(i, items);
                                }
                            }
                        }
                    }
                }
                state_val.items.push(BedGraphSectionItem {
                    start: current_val.start,
                    end: current_val.end,
                    val: current_val.value,
                });

                match &mut summary {
                    None => {
                        summary = Some(Summary {
                            bases_covered: (current_val.end - current_val.start) as u64,
                            min_val: current_val.value as f64,
                            max_val: current_val.value as f64,
                            sum: (current_val.end - current_val.start) as f64 * current_val.value as f64,
                            sum_squares: (current_val.end - current_val.start) as f64 * (current_val.value * current_val.value) as f64,
                        })
                    },
                    Some(summary) => {
                        summary.bases_covered += (current_val.end - current_val.start) as u64;
                        summary.min_val = summary.min_val.min(current_val.value as f64);
                        summary.max_val = summary.max_val.max(current_val.value as f64);
                        summary.sum += (current_val.end - current_val.start) as f64 * current_val.value as f64;
                        summary.sum_squares += (current_val.end - current_val.start) as f64 * (current_val.value * current_val.value) as f64;
                    }
                }
            }

            // We have previously had data, need to write final section (if needed)
            if let Some(state_val) = state {
                let lastchrom = state_val.chrom;
                let chromId = chrom_ids.get_id(lastchrom.clone());
                if state_val.items.len() > 0 {
                    write_section_items(state_val.items, chromId, &mut ftx);
                }

                for (i, mut zoom_item) in state_val.zoom_items.into_iter().enumerate() {
                    match &zoom_item.live_info {
                        None => (),
                        Some(zoom2) => {
                            assert!(lastchrom == zoom2.chrom);
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
                            //println!("Zoom file size ({:?}): {:?}", zoom.0, zoom.1.seek(SeekFrom::End(0))?)
                        }
                    }
                    if zoom_item.records.len() > 0 {
                        write_zoom_items(i, zoom_item.records);
                    }
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
            Ok((chrom_ids, summary_complete))
        };

        let f = read_file();

        let process_zooms = async move || -> std::io::Result<Vec<TempZoomInfo>> {
            let mut current_offsets: Vec<u64> = (0..num_zooms).map(|_| 0).collect();
            while let Some((idx, future)) = await!(frxzooms.next()) {
                let section_raw = await!(future);
                let section = section_raw?;
                let size = section.data.len() as u64;
                let file = &mut zooms[idx as usize].1;
                let current_offset = &mut current_offsets[idx as usize];
                file.write_all(&section.data)?;
                zoom_sections_out[idx as usize].push(Section {
                    chrom: section.chrom,
                    start: section.start,
                    end: section.end,
                    offset: *current_offset,
                    size,
                });
                *current_offset += size;
            }
            let zooms_info = zooms.into_iter().zip(zoom_sections_out.into_iter())
                .map(|(z, sections)|
                    (z.0, z.1.into_inner().expect("Can't get inner file"), sections)
                    ).collect();
            Ok(zooms_info)
        };

        let zooms_future = process_zooms();

        let (zooms_remote, zooms_handle) = zooms_future.remote_handle();
        std::thread::spawn(move || {
            futures::executor::block_on(zooms_remote);
        });

        let (sections_remote, sections_handle) = sections_future.remote_handle();
        std::thread::spawn(move || {
            futures::executor::block_on(sections_remote);
        });

        Ok((f.map(|ab| ab.unwrap()), frxsections, zooms_handle.map(|z| z.unwrap()), sections_handle.map(|f| f.unwrap())))
    }

    fn write_zooms(&mut self, zooms: Vec<TempZoomInfo>, zoom_entries: &mut Vec<ZoomHeader>) -> std::io::Result<()> {
        // TODO: param
        const DO_COMPRESS: bool = true;
        const ITEMS_PER_SLOT: u32 = 1024;

        for mut zoom in zooms {
            let zoom_data_offset = {
                let mut file = self.get_buf_writer();
                file.seek(SeekFrom::Current(0))?
            };
            let (sections_stream, zoom_index_offset) = {
                let mut file = self.get_buf_writer();
                zoom.1.seek(SeekFrom::Start(0))?;
                let sections_iter = zoom.2.iter().map(|section| {
                    let chrom = section.chrom;
                    let start = section.start;
                    let end = section.end;
                    Section {
                        offset: zoom_data_offset + section.offset,
                        size: section.size,
                        chrom: chrom,
                        start: start,
                        end: end,
                    }
                });
                let sections_stream = futures::stream::iter(sections_iter);

                let mut buf_reader = std::io::BufReader::new(zoom.1);
                std::io::copy(&mut buf_reader, &mut file)?;
                let zoom_index_offset = file.seek(SeekFrom::Current(0))?;

                (sections_stream, zoom_index_offset)
            };
            let (nodes, levels, total_sections) = futures::executor::block_on(BigWig::get_rtreeindex(sections_stream));
            self.write_rtreeindex(nodes, levels, total_sections)?;

            zoom_entries.push(ZoomHeader {
                reduction_level: zoom.0,
                data_offset: zoom_data_offset,
                index_offset: zoom_index_offset,
            });
        }

        Ok(())
    }

    async fn get_rtreeindex<S>(sections_stream: S) -> (RTreeNodeList<RTreeNode>, usize, u64) where S : Stream<Item=Section> + Unpin {
        const BLOCK_SIZE: u32 = 256;

        let mut total_sections = 0;
        let mut current_nodes: Box<Stream<Item=RTreeNode> + Unpin> = Box::new(sections_stream.map(|s| RTreeNode {
            start_chrom_idx: s.chrom,
            start_base: s.start,
            end_chrom_idx: s.chrom,
            end_base: s.end,
            kind: RTreeNodeType::Leaf {
                offset: s.offset,
                size: s.size,
            },
        }));
        let mut levels = 0;
        let nodes: RTreeNodeList<RTreeNode> = loop {
            let mut start_chrom_idx = 0;
            let mut start_base = 0;
            let mut end_chrom_idx = 0;
            let mut end_base = 0;
            let mut next_nodes: Vec<RTreeNode> = vec![];
            let mut current_group: Vec<RTreeNode> = vec![];
            loop {
                let next_node = await!(current_nodes.next());
                match next_node {
                    None => {
                        //println!("Remaining nodes at complete: {}", current_group.len());
                        if next_nodes.len() == 0 {
                            next_nodes = current_group;
                        } else {
                            levels += 1;
                            next_nodes.push(RTreeNode{
                                start_chrom_idx,
                                start_base,
                                end_chrom_idx,
                                end_base,
                                kind: RTreeNodeType::NonLeaf {
                                    children: RTreeNodeList::<RTreeNode> {
                                        nodes: current_group
                                    }
                                },
                            });
                        }
                        break
                    },
                    Some(node) => {
                        if levels == 0 {
                            total_sections += 1;
                        }
                        if current_group.len() == 0 {
                            start_chrom_idx = node.start_chrom_idx;
                            start_base = node.start_base;
                            end_chrom_idx = node.end_chrom_idx;
                            end_base = node.end_base;
                        }
                        if end_chrom_idx == node.end_chrom_idx {
                            end_base = std::cmp::max(end_base, node.end_base);
                        } else {
                            end_base = node.end_base
                        }
                        end_chrom_idx = std::cmp::max(end_chrom_idx, node.end_chrom_idx);
                        current_group.push(node);
                        if current_group.len() >= BLOCK_SIZE as usize {
                            next_nodes.push(RTreeNode{
                                start_chrom_idx,
                                start_base,
                                end_chrom_idx,
                                end_base,
                                kind: RTreeNodeType::NonLeaf {
                                    children: RTreeNodeList::<RTreeNode> {
                                        nodes: current_group
                                    }
                                },
                            });
                            current_group = vec![];
                        }
                    }
                }
            }            

            if next_nodes.len() < BLOCK_SIZE as usize {
                break RTreeNodeList::<RTreeNode> {
                    nodes: next_nodes
                }
            }

            current_nodes = Box::new(stream::iter(next_nodes.into_iter()));
        };
        //println!("Total sections: {:?}", total_sections);
        //println!("Nodes ({:?}): {:?}", nodes.nodes.len(), nodes);
        //println!("Levels: {:?}", levels);
        (nodes, levels, total_sections)
    }

    fn write_rtreeindex(&mut self, nodes: RTreeNodeList<RTreeNode>, levels: usize, section_count: u64) -> std::io::Result<()> {
        const BLOCK_SIZE: u32 = 256;

        const NODEHEADER_SIZE: u64 = 1 + 1 + 2;
        const NON_LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8;
        const LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8 + 8;

        let mut file = self.get_buf_writer();

        let mut index_offsets: Vec<u64> = vec![0u64; levels as usize];

        fn calculate_offsets(mut index_offsets: &mut Vec<u64>, trees: &RTreeNodeList<RTreeNode>, level: usize) -> std::io::Result<()> {
            if level == 0 {
                return Ok(())
            }
            let isleaf: bool = {
                if trees.nodes.len() == 0 {
                    false
                } else {
                    match trees.nodes[0].kind {
                        RTreeNodeType::Leaf { .. } => true,
                        RTreeNodeType::NonLeaf { .. } => false,
                    }
                }
            };
            index_offsets[level - 1] += NODEHEADER_SIZE;
            for tree in trees.nodes.iter() {
                match &tree.kind {
                    RTreeNodeType::Leaf { .. } => panic!("Only calculating offsets/sizes for indices (level > 0)"),
                    RTreeNodeType::NonLeaf { children, .. } => {
                        debug_assert!(level != 0, "Non Leaf node found at level 0");
                        debug_assert!(!isleaf, "Mixed node types at level {}", level);

                        index_offsets[level - 1] += NON_LEAFNODE_SIZE;

                        calculate_offsets(&mut index_offsets, &children, level - 1)?;
                    },
                }
            }
            Ok(())
        }

        if levels > 0 {
            calculate_offsets(&mut index_offsets, &nodes, levels)?;
        }
        //println!("index Offsets: {:?}", index_offsets);


        const NON_LEAFNODE_FULL_BLOCK_SIZE: u64 = NODEHEADER_SIZE + NON_LEAFNODE_SIZE * BLOCK_SIZE as u64;
        const LEAFNODE_FULL_BLOCK_SIZE: u64 = NODEHEADER_SIZE + LEAFNODE_SIZE * BLOCK_SIZE as u64;
        fn write_tree(mut file: &mut BufWriter<&File>, index_offsets: &Vec<u64>, trees: &RTreeNodeList<RTreeNode>, curr_level: usize, dest_level: usize, childnode_offset: u64) -> std::io::Result<u64> {
            assert!(curr_level >= dest_level);
            let mut total_size = 0;
            if curr_level != dest_level {
                let mut next_offset_offset = 0;
                for tree in trees.nodes.iter() {
                    match &tree.kind {
                        RTreeNodeType::Leaf { .. } => panic!("Leaf node found at level {}", curr_level),
                        RTreeNodeType::NonLeaf { children, .. } => {
                            next_offset_offset += write_tree(&mut file, &index_offsets, &children, curr_level - 1, dest_level, childnode_offset + next_offset_offset)?;
                        },
                    }
                }
                total_size += next_offset_offset;
                return Ok(total_size)
            }
            let isleaf = if let RTreeNodeType::Leaf { .. } = trees.nodes[0].kind {
                1
            } else {
                0
            };

            //println!("Writing {}. Isleaf: {} At: {}", trees.nodes.len(), isleaf, file.seek(SeekFrom::Current(0))?);
            file.write_u8(isleaf)?;
            file.write_u8(0)?;
            file.write_u16::<NativeEndian>(trees.nodes.len() as u16)?;
            total_size += 4;
            for (idx, node) in trees.nodes.iter().enumerate() {
                file.write_u32::<NativeEndian>(node.start_chrom_idx)?;
                file.write_u32::<NativeEndian>(node.start_base)?;
                file.write_u32::<NativeEndian>(node.end_chrom_idx)?;
                file.write_u32::<NativeEndian>(node.end_base)?;
                total_size += 16;
                match &node.kind {
                    RTreeNodeType::Leaf { offset, size } => {
                        file.write_u64::<NativeEndian>(*offset)?;
                        file.write_u64::<NativeEndian>(*size)?;
                        total_size += 16;
                    },
                    RTreeNodeType::NonLeaf { .. } => {
                        debug_assert!(curr_level != 0);
                        let full_size = if (curr_level - 1) > 0 {
                            NON_LEAFNODE_FULL_BLOCK_SIZE
                        } else {
                            LEAFNODE_FULL_BLOCK_SIZE
                        };
                        let child_offset: u64 = childnode_offset + idx as u64 * full_size;
                        //println!("Child node offset: {}; Added: {}", child_offset, idx as u64 * full_size);
                        file.write_u64::<NativeEndian>(child_offset)?;
                        total_size += 8;
                    },
                }
            }
            Ok(total_size)
        }

        // TODO remove
        const ITEMS_PER_SLOT: u32 = 1024;

        let end_of_data = file.seek(SeekFrom::Current(0))?;
        // TODO: handle case with no data
        {
            //println!("cirTree header (write):\n bs: {:?}\n ic: {:?}\n sci: {:?}\n sb: {:?}\n eci: {:?}\n eb: {:?}\n efo: {:?}\n ips: {:?}\n r: {:?}", BLOCK_SIZE, section_count, nodes.nodes[0].start_chrom_idx, nodes.nodes[0].start_base, nodes.nodes[nodes.nodes.len() - 1].end_chrom_idx, nodes.nodes[nodes.nodes.len() - 1].end_base, end_of_data, ITEMS_PER_SLOT, 0);
            file.write_u32::<NativeEndian>(CIR_TREE_MAGIC)?;
            file.write_u32::<NativeEndian>(BLOCK_SIZE)?;
            file.write_u64::<NativeEndian>(section_count)?;
            file.write_u32::<NativeEndian>(nodes.nodes[0].start_chrom_idx)?;
            file.write_u32::<NativeEndian>(nodes.nodes[0].start_base)?;
            file.write_u32::<NativeEndian>(nodes.nodes[nodes.nodes.len() - 1].end_chrom_idx)?;
            file.write_u32::<NativeEndian>(nodes.nodes[nodes.nodes.len() - 1].end_base)?;
            file.write_u64::<NativeEndian>(end_of_data)?;
            file.write_u32::<NativeEndian>(ITEMS_PER_SLOT)?;
            file.write_u32::<NativeEndian>(0)?;
        }

        let mut current_offset = file.seek(SeekFrom::Current(0))?;
        //println!("Levels: {:?}", levels);
        //println!("Start of index: {}", current_offset);
        for level in (0..=levels).rev() {
            if level > 0 {
                current_offset += index_offsets[level - 1];
            }
            write_tree(&mut file, &index_offsets, &nodes, levels, level, current_offset)?;
            //println!("End of index level {}: {}", level, file.seek(SeekFrom::Current(0))?);
        }

        Ok(())
    }
}
