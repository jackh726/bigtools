#![allow(non_snake_case)]
#![allow(dead_code)]


extern crate flate2;
extern crate tempfile;

use std::vec::Vec;

use std::io::{Seek, SeekFrom};
use std::io::prelude::*;
use std::fs::File;
use std::io::BufWriter;

use byteorder::{LittleEndian, BigEndian, NativeEndian, ReadBytesExt, WriteBytesExt};

use flate2::Compression;
use flate2::write::ZlibEncoder;
use flate2::read::ZlibDecoder;

const BIGWIG_MAGIC_LTH: u32 = 0x888FFC26;
const BIGWIG_MAGIC_HTL: u32 = 0x26FC8F88;
const BIGBED_MAGIC_LTH: u32 = 0x8789F2EB;
const BIGBED_MAGIC_HTL: u32 = 0xEBF28987;

const CIR_TREE_MAGIC: u32 = 0x2468ACE0;
const CHROM_TREE_MAGIC: u32 = 0x78CA8C91;

#[derive(Debug)]
struct BBIHeader {
    isBigEndian: bool,

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

type TempZoomInfo = (u32 /* resolution */, File /* Temp file that contains data */, u32 /* Number of zoom records */);

struct BedGraphSectionItem {
    start: u32,
    end: u32,
    val: f32,
}

pub struct BigWig {
    pub path: String,
    fp: File,
    info: Option<Box<BigWigInfo>>,
}

impl BigWig {
    pub fn test_read_zoom(&mut self, chrom_name: &str, start: u32, end: u32) -> std::io::Result<()> {
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
        {
            let file = &mut self.fp;
            let info = self.info.as_ref().unwrap();
            file.seek(SeekFrom::Start(info.zoom_headers[0].index_offset))?;
        }
        let blocks = self.get_overlapping_blocks(chrom_ix, start, end)?;

        for block in blocks {
            let info = self.info.as_ref().unwrap();
            let bigendian = info.header.isBigEndian;
            let file = &mut self.fp;
            file.seek(SeekFrom::Start(block.offset))?;

            let mut raw_data = vec![0u8; block.size as usize];
            file.read_exact(&mut raw_data)?;
            let data = if info.header.uncompress_buf_size > 0 {
                let mut uncompressed_block_data = vec![0u8; info.header.uncompress_buf_size as usize];
                let mut d = ZlibDecoder::new(&raw_data[..]);
                d.read(&mut uncompressed_block_data)?;
                uncompressed_block_data
            } else {
                raw_data
            };
            let itemcount = data.len() / (4 * 8);
            assert!(data.len() % (4 * 8) == 0);
            let mut data_mut = &data[..];
            for _ in 0..itemcount {
                read_u32(&mut data_mut, bigendian)?;
                read_u32(&mut data_mut, bigendian)?;
                read_u32(&mut data_mut, bigendian)?;
                read_u32(&mut data_mut, bigendian)?;
                read_f32(&mut data_mut, bigendian)?;
                read_f32(&mut data_mut, bigendian)?;
                read_f32(&mut data_mut, bigendian)?;
                read_f32(&mut data_mut, bigendian)?;
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

        let header = {
            let mut file = &mut self.fp;

            let magic = read_u32(&mut file, false)?;
            println!("Magic {:x?}: ", magic);
            let bigendian = match magic {
                BIGWIG_MAGIC_HTL => true,
                BIGWIG_MAGIC_LTH => false,
                _ => return Err(std::io::Error::new(std::io::ErrorKind::Other, "File not a big wig"))
            };

            let version = read_u16(&mut file, bigendian)?;
            let zoom_levels = read_u16(&mut file, bigendian)?;
            let chromosome_tree_offset = read_u64(&mut file, bigendian)?;
            let full_data_offset = read_u64(&mut file, bigendian)?;
            let full_index_offset = read_u64(&mut file, bigendian)?;
            let field_count = read_u16(&mut file, bigendian)?;
            let defined_field_count = read_u16(&mut file, bigendian)?;
            let auto_sql_offset = read_u64(&mut file, bigendian)?;
            let total_summary_offset = read_u64(&mut file, bigendian)?;
            let uncompress_buf_size = read_u32(&mut file, bigendian)?;
            let reserved = read_u64(&mut file, bigendian)?;

            BBIHeader {
                isBigEndian: bigendian,
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
            }
        };
        println!("Header: {:?}", header);

        let zoom_headers = self.read_zoom_headers(&header)?;

        let chrom_info = {
            let mut file = &mut self.fp;

            file.seek(SeekFrom::Start(header.chromosome_tree_offset))?;
            let magic = read_u32(&mut file, header.isBigEndian)?;
            let _block_size = read_u32(&mut file, header.isBigEndian)?;
            let key_size = read_u32(&mut file, header.isBigEndian)?;
            let val_size = read_u32(&mut file, header.isBigEndian)?;
            let item_count = read_u64(&mut file, header.isBigEndian)?;
            let _reserved = read_u64(&mut file, header.isBigEndian)?;
            if magic != CHROM_TREE_MAGIC {
                return Err(std::io::Error::new(std::io::ErrorKind::Other, "Invalid file format: CHROM_TREE_MAGIC does not match."))
            }
            //println!("{:x?} {:?} {:?} {:?} {:?} {:?}", magic, _block_size, key_size, val_size, item_count, _reserved);
            assert_eq!(val_size, 8u32); 

            let mut chroms = Vec::with_capacity(item_count as usize);
            BigWig::read_chrom_tree_block(&mut file, header.isBigEndian, &mut chroms, key_size)?;
            chroms
        };

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

    fn read_zoom_headers(&mut self, header: &BBIHeader) -> std::io::Result<Vec<ZoomHeader>> {
        let mut file = &mut self.fp;

        let mut zoom_headers = vec![];
        for _ in 1..header.zoom_levels {
            let reduction_level = read_u32(&mut file, header.isBigEndian)?;
            let _reserved = read_u32(&mut file, header.isBigEndian)?;
            let data_offset = read_u64(&mut file, header.isBigEndian)?;
            let index_offset = read_u64(&mut file, header.isBigEndian)?;

            println!("Zoom header: reductionLevel: {:?} Reserved: {:?} Data offset: {:?} Index offset: {:?}", reduction_level, _reserved, data_offset, index_offset);

            zoom_headers.push(ZoomHeader {
                reduction_level,
                data_offset,
                index_offset,
            });
        }

        Ok(zoom_headers)
    }

    fn read_chrom_tree_block(mut f: &mut File, bigendian: bool, chroms: &mut Vec<ChromInfo>, key_size: u32) -> std::io::Result<()> {
        let isleaf = read_u8(&mut f, bigendian)?;
        let _reserved = read_u8(&mut f, bigendian)?;
        let count = read_u16(&mut f, bigendian)?;

        if isleaf == 1 {
            for _ in 0..count {
                let mut key_bytes = vec![0u8; key_size as usize];
                f.read_exact(&mut key_bytes)?;
                let key_string = match String::from_utf8(key_bytes) {
                    Ok(s) => s.trim_matches(char::from(0)).to_owned(),
                    Err(_) => return Err(std::io::Error::new(std::io::ErrorKind::Other, "Invalid file format: Invalid utf-8 string.")),
                };
                let chrom_id = read_u32(&mut f, bigendian)?;
                let chrom_size = read_u32(&mut f, bigendian)?;
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
                let child_offset = read_u64(&mut f, bigendian)?;
                current_position = f.seek(SeekFrom::Current(0))?;
                f.seek(SeekFrom::Start(child_offset))?;
                BigWig::read_chrom_tree_block(f, bigendian, chroms, key_size)?;
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

    fn search_overlapping_blocks(&mut self, chrom_ix: u32, start: u32, end: u32, mut blocks: &mut Vec<Block>) -> std::io::Result<()> {
        println!("Searching for overlapping blocks at {:?}. Searching {:?}:{:?}-{:?}", self.current_file_offset()?, chrom_ix, start, end);

        let bigendian = {
            let info = self.info.as_ref().unwrap();
            (&info.header).isBigEndian
        };
        let (isleaf, count): (u8, u16) = {
            let mut file = self.get_buf_reader();
            let isleaf = read_u8(&mut file, bigendian)?;
            let _reserved = read_u8(&mut file, bigendian)?;
            let count = read_u16(&mut file, bigendian)?;
            println!("Index: {:?} {:?} {:?}", isleaf, _reserved, count);
            (isleaf, count)
        };
        if isleaf == 1 {
            for _ in 0..count {
                let mut file = self.get_buf_reader();
                let start_chrom_ix = read_u32(&mut file, bigendian)?;
                let start_base = read_u32(&mut file, bigendian)?;
                let end_chrom_ix = read_u32(&mut file, bigendian)?;
                let end_base = read_u32(&mut file, bigendian)?;
                let data_offset = read_u64(&mut file, bigendian)?;
                let data_size = read_u64(&mut file, bigendian)?;
                if !BigWig::overlaps(chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base) {
                    continue;
                }
                println!("Overlaps (leaf): {:?}:{:?}-{:?} with {:?}:{:?}-{:?}:{:?} {:?} {:?}", chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset, data_size);
                blocks.push(Block {
                    offset: data_offset,
                    size: data_size,
                })
            }
        } else {
            let mut childblocks: Vec<u64> = vec![];
            {
                let mut file = self.get_buf_reader();
                for _ in 0..count {
                    let start_chrom_ix = read_u32(&mut file, bigendian)?;
                    let start_base = read_u32(&mut file, bigendian)?;
                    let end_chrom_ix = read_u32(&mut file, bigendian)?;
                    let end_base = read_u32(&mut file, bigendian)?;
                    let data_offset = read_u64(&mut file, bigendian)?;
                    if !BigWig::overlaps(chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base) {
                        continue;
                    }
                    println!("Overlaps (non-leaf): {:?}:{:?}-{:?} with {:?}:{:?}-{:?}:{:?} {:?}", chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset);
                    childblocks.push(data_offset);
                }
            }
            for childblock in childblocks {
                println!("Seeking to {:?}", childblock);
                (&mut self.fp).seek(SeekFrom::Start(childblock))?;
                self.search_overlapping_blocks(chrom_ix, start, end, &mut blocks)?;
            }
        }
        return Ok(());
    }

    fn get_overlapping_blocks(&mut self, chrom_ix: u32, start: u32, end: u32) -> std::io::Result<Vec<Block>> {
        let bigendian = {
            let info = self.info.as_ref().unwrap();
            (&info.header).isBigEndian
        };
        let blocks = {
            let mut file = self.get_buf_reader();
            let magic = read_u32(&mut file, bigendian)?;
            if magic != CIR_TREE_MAGIC {
                return Err(std::io::Error::new(std::io::ErrorKind::Other, "Invalid file format: CIR_TREE_MAGIC does not match."));
            }
            let _blocksize = read_u32(&mut file, bigendian)?;
            let _item_count = read_u64(&mut file, bigendian)?;
            let _start_chrom_idx = read_u32(&mut file, bigendian)?;
            let _start_base = read_u32(&mut file, bigendian)?;
            let _end_chrom_idx = read_u32(&mut file, bigendian)?;
            let _end_base = read_u32(&mut file, bigendian)?;
            let _end_file_offset = read_u64(&mut file, bigendian)?;
            let _item_per_slot = read_u32(&mut file, bigendian)?;
            let _reserved = read_u32(&mut file, bigendian)?;

            // TODO: could do some optimization here to check if our interval overlaps with any data

            println!("cirTree header:\n bs: {:?}\n ic: {:?}\n sci: {:?}\n sb: {:?}\n eci: {:?}\n eb: {:?}\n efo: {:?}\n ips: {:?}\n r: {:?}", _blocksize, _item_count, _start_chrom_idx, _start_base, _end_chrom_idx, _end_base, _end_file_offset, _item_per_slot, _reserved);
            let mut blocks: Vec<Block> = vec![];
            self.search_overlapping_blocks(chrom_ix, start, end, &mut blocks)?;
            println!("overlapping_blocks: {:?}", blocks);
            blocks
        };
        Ok(blocks)
    }

    pub fn get_interval(&mut self, chrom_name: &str, start: u32, end: u32) -> std::io::Result<Vec<Value>> {
        self.ensure_info().or(Err(std::io::Error::new(std::io::ErrorKind::Other, "Must first call read_info")))?;
        let bigendian = {
            let info = self.info.as_ref().unwrap();
            (&info.header).isBigEndian
        };
        let blocks = {
            {
                let info = self.info.as_ref().unwrap();
                let file = &mut self.fp;
                file.seek(SeekFrom::Start((&info.header).full_index_offset))?;
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
            self.get_overlapping_blocks(chrom_ix, start, end)?
        };

        let mut values: Vec<Value> = Vec::new();
        let uncompress_buf_size: usize = {
            let info = self.info.as_ref().unwrap();
            info.header.uncompress_buf_size as usize
        };
        let mut file = self.get_buf_reader();
        for block in blocks {
            // TODO: Could minimize this by chunking block reads
            let block_data: Vec<u8> = {
                file.seek(SeekFrom::Start(block.offset))?;
                let mut raw_data = vec![0u8; block.size as usize];
                file.read_exact(&mut raw_data)?;
                let data = if uncompress_buf_size > 0 {
                    let mut uncompressed_block_data = vec![0u8; uncompress_buf_size];
                    let mut d = ZlibDecoder::new(&raw_data[..]);
                    d.read(&mut uncompressed_block_data)?;
                    uncompressed_block_data
                } else {
                    raw_data
                };
                data
            };
            let mut block_data_mut = &block_data[..];
            let _chrom_id = read_u32(&mut block_data_mut, bigendian)?;
            let chrom_start = read_u32(&mut block_data_mut, bigendian)?;
            let _chrom_end = read_u32(&mut block_data_mut, bigendian)?;
            let item_step = read_u32(&mut block_data_mut, bigendian)?;
            let item_span = read_u32(&mut block_data_mut, bigendian)?;
            let section_type = read_u8(&mut block_data_mut, bigendian)?;
            let _reserved = read_u8(&mut block_data_mut, bigendian)?;
            let item_count = read_u16(&mut block_data_mut, bigendian)?;

            let mut start = chrom_start;
            for _ in 0..item_count {
                match section_type {
                    1 => {
                        // bedgraph
                        let chrom_start = read_u32(&mut block_data_mut, bigendian)?;
                        let chrom_end = read_u32(&mut block_data_mut, bigendian)?;
                        let value = read_f32(&mut block_data_mut, bigendian)?;
                        values.push(Value {
                            start: std::cmp::max(chrom_start, start),
                            end: std::cmp::min(chrom_end, end),
                            value,
                        });
                    },
                    2 => {
                        // variable step
                        let chrom_start = read_u32(&mut block_data_mut, bigendian)?;
                        let chrom_end = chrom_start + item_span;
                        let value = read_f32(&mut block_data_mut, bigendian)?;
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
                        let value = read_f32(&mut block_data_mut, bigendian)?;
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

    pub fn write<'a, V>(&mut self, chrom_sizes: std::collections::HashMap<&str, u32>, mut vals: V) -> std::io::Result<()> where V : std::iter::Iterator<Item=ValueWithChrom> {
        self.write_blank_headers()?;

        let total_summary_offset = self.current_file_offset()?;

        self.write_blank_summary()?;

        let mut chrom_ids = std::collections::HashMap::new();

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
        let (sections, summary, zooms) = self.write_vals(&mut vals, &mut chrom_ids)?;
        let total_sections = sections.len() as u32;
        let data_size = self.current_file_offset()? - pre_data;
        println!("Data size: {:?}", data_size);
        println!("Sections: {:?}", sections.len());
        println!("Summary: {:?}", summary);
        println!("Zooms: {:?}", zooms);

        // Since the chrom tree is read before the index, we put this before the full data index
        // Therefore, there is a higher likelihood that the udc file will only need one read for chrom tree + full data index
        // Putting the chrom tree before the data also has a higher likelihood of being included with the beginning headers,
        //  but requires us to know all the data ahead of time (when writing)
        let chrom_index_start = self.current_file_offset()?;
        self.write_chrom_tree(chrom_sizes, &chrom_ids)?;

        let index_start = self.current_file_offset()?;
        self.write_rtreeindex(sections)?;

        const DO_COMPRESS: bool = true; // TODO: param
        const ITEMS_PER_SLOT: u32 = 1024; // TODO: param

        let mut zoom_entries: Vec<ZoomHeader> = vec![];
        self.write_zooms(zooms, &mut zoom_entries)?;

        println!("Zoom entries: {:?}", zoom_entries);
        let num_zooms = zoom_entries.len() as u16;

        // We *could* actually check the the real max size, but let's just assume at it's as large as the largest possible value
        // In most cases, I think this is the true max size (unless there is only one section and its less than ITEMS_PER_SLOT in size)
        let uncompress_buf_size = if DO_COMPRESS {
            ITEMS_PER_SLOT * (1 + 1 + 2 + 4 + 4 + 4 + 4 + 8 + 8)
        } else {
            0
        };
        {
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
        }
        {
            let mut file = self.get_buf_writer();
            for zoom_entry in zoom_entries {
                file.write_u32::<NativeEndian>(zoom_entry.reduction_level)?;
                file.write_u32::<NativeEndian>(0)?;
                file.write_u64::<NativeEndian>(zoom_entry.data_offset)?;
                file.write_u64::<NativeEndian>(zoom_entry.index_offset)?;
            }
        }
        {
            let mut file = self.get_buf_writer();
            assert!(file.seek(SeekFrom::Current(0))? == total_summary_offset);
            file.write_u64::<NativeEndian>(summary.bases_covered)?;
            file.write_f64::<NativeEndian>(summary.min_val)?;
            file.write_f64::<NativeEndian>(summary.max_val)?;
            file.write_f64::<NativeEndian>(summary.sum)?;
            file.write_f64::<NativeEndian>(summary.sum_squares)?;
        }
        {
            let mut file = self.get_buf_writer();
            file.write_u32::<NativeEndian>(total_sections)?;

            file.seek(SeekFrom::End(0))?;
            file.write_u32::<NativeEndian>(BIGWIG_MAGIC_LTH)?;
        }

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

    fn write_chrom_tree(&mut self, chrom_sizes: std::collections::HashMap<&str, u32>, chrom_ids: &std::collections::HashMap<String, u32>) -> std::io::Result<()> {
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
        println!("Max bytes: {:?}", max_bytes);
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
            println!("Key bytes for {:?}: {:x?}", chrom, key_bytes);
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

    fn write_section(file: &mut std::io::BufWriter<&File>, items_in_section: Vec<BedGraphSectionItem>, chromId: u32, current_offset: u64) -> std::io::Result<(Section, u64)> {
        let mut bytes: Vec<u8> = vec![];
        
        //let current_offset = file.seek(SeekFrom::Current(0))?;
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
        let out_length = out_bytes.len();
        file.write_all(&out_bytes)?;
        //let end_offset = file.seek(SeekFrom::Current(0))?;
        let end_offset = current_offset + out_length as u64;


        Ok((Section {
            offset: current_offset,
            size: end_offset - current_offset,
            chrom: chromId,
            start,
            end,
        }, end_offset))
    }

    fn write_vals<'a, V>(&mut self, vals_iter: &mut V, chrom_ids: &mut std::collections::HashMap<String, u32>) -> std::io::Result<(Vec<Section>, Summary, Vec<TempZoomInfo>)> where V : std::iter::Iterator<Item=ValueWithChrom> {
        let ITEMS_PER_SLOT: u16 = 1024;

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

        type ZoomInfo<'a> = (u32, BufWriter<File>, Option<LiveZoomInfo>, u32);

        let zoom_sizes: Vec<u32> = vec![10, 40, 160, 640, 2560, 10240, 40960, 163840, 655360, 2621440];
        //let zoom_sizes: Vec<u32> = vec![];
        let mut zooms: Vec<ZoomInfo> = vec![];
        for size in zoom_sizes {
            let temp = tempfile::tempfile()?;
            zooms.push((size, std::io::BufWriter::new(temp), None, 0));
        }
        
        let mut next_chrom_id: u32 = 0;
        let mut sections = Vec::new();
        let mut summary: Option<Summary> = None;
        let mut last_chrom: Option<String> = None;
        let mut last_start: u32 = 0;
        let mut items_in_section: Vec<BedGraphSectionItem> = Vec::with_capacity(ITEMS_PER_SLOT as usize);

        let write_zoom = |file: &mut Write, zoom_info: &Option<LiveZoomInfo>, chromId: u32| -> std::io::Result<()> {
            match &zoom_info {
                None => (),
                Some(live) => {
                    file.write_u32::<NativeEndian>(chromId)?;
                    file.write_u32::<NativeEndian>(live.start)?;
                    file.write_u32::<NativeEndian>(live.end)?;
                    file.write_u32::<NativeEndian>(live.valid_count)?;
                    file.write_f32::<NativeEndian>(live.min_value)?;
                    file.write_f32::<NativeEndian>(live.max_value)?;
                    file.write_f32::<NativeEndian>(live.sum)?;
                    file.write_f32::<NativeEndian>(live.sum_squares)?;
                }
            }
            Ok(())
        };

        let file = &mut self.get_buf_writer();
        let mut current_offset = file.seek(SeekFrom::Current(0))?;
        for current_val in vals_iter {
            if let Some(lastchrom) = &last_chrom {
                assert!(lastchrom <= &current_val.chrom);
                if lastchrom == &current_val.chrom {
                    assert!(last_start <= current_val.start);
                }
            }
            assert!(current_val.start <= current_val.end);

            if items_in_section.len() >= ITEMS_PER_SLOT as usize || (last_chrom.is_some() && &current_val.chrom != last_chrom.as_ref().unwrap()) {
                let chrom = last_chrom.unwrap();
                let chromId: u32 = *chrom_ids.entry(chrom.clone()).or_insert(next_chrom_id);
                if chromId == next_chrom_id {
                    next_chrom_id += 1;
                }

                let (section, next_offset) = BigWig::write_section(file, items_in_section, chromId, current_offset)?;
                current_offset = next_offset;
                sections.push(section);
                items_in_section = vec![];

                if current_val.chrom != chrom {
                    for mut zoom in &mut zooms {
                        write_zoom(&mut zoom.1, &zoom.2, chromId)?;
                        zoom.3 += 1;
                        zoom.2 = None;
                    }
                }
            }
            for mut zoom in &mut zooms {
                let mut add_start = current_val.start;
                loop {
                    match &mut zoom.2 {
                        None => {
                            zoom.2 = Some(LiveZoomInfo {
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
                            if zoom2.end >= current_val.end {
                                break;
                            }
                            let next_end = zoom2.start + zoom.0;
                            if next_end < current_val.start {
                                // The last zoom entry ended before this value begins, need to write
                                let chromId: u32 = *chrom_ids.entry(zoom2.chrom.to_string()).or_insert(next_chrom_id);
                                if chromId == next_chrom_id {
                                    next_chrom_id += 1;
                                }
                                write_zoom(&mut zoom.1, &zoom.2, chromId)?;
                                zoom.3 += 1;
                                zoom.2 = None;
                                continue;
                            }
                            if next_end >= current_val.end {
                                // The current value isn't enough to finish this zoom entry
                                let added_bases = current_val.end - add_start;                                
                                zoom2.end = current_val.end;
                                zoom2.valid_count += added_bases;
                                zoom2.min_value = zoom2.min_value.min(current_val.value);
                                zoom2.max_value = zoom2.max_value.max(current_val.value);
                                zoom2.sum += added_bases as f32 * current_val.value;
                                zoom2.sum_squares += added_bases as f32 * current_val.value * current_val.value;
                                if next_end == current_val.end {
                                    let chromId: u32 = *chrom_ids.entry(zoom2.chrom.to_string()).or_insert(next_chrom_id);
                                    if chromId == next_chrom_id {
                                        next_chrom_id += 1;
                                    }
                                    write_zoom(&mut zoom.1, &zoom.2, chromId)?;
                                    zoom.3 += 1;
                                    zoom.2 = None;
                                    break;
                                }
                            } else {
                                // The current value will finish the zoom entry
                                let added_bases = next_end - add_start;
                                zoom2.end = next_end;
                                zoom2.valid_count += added_bases;
                                zoom2.min_value = zoom2.min_value.min(current_val.value);
                                zoom2.max_value = zoom2.max_value.max(current_val.value);
                                zoom2.sum += added_bases as f32 * current_val.value;
                                zoom2.sum_squares += added_bases as f32 * current_val.value * current_val.value;
                                let chromId: u32 = *chrom_ids.entry(zoom2.chrom.to_string()).or_insert(next_chrom_id);
                                if chromId == next_chrom_id {
                                    next_chrom_id += 1;
                                }
                                write_zoom(&mut zoom.1, &zoom.2, chromId)?;
                                zoom.3 += 1;
                                zoom.2 = None;
                                add_start = next_end;
                            }
                        }
                    }

                }
            }
            items_in_section.push(BedGraphSectionItem {
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

            last_start = current_val.start;
            last_chrom = Some(current_val.chrom);
        }
        if let Some(lastchrom) = last_chrom {
            let chromId: u32 = *chrom_ids.entry(lastchrom.to_string()).or_insert(next_chrom_id);
            //if chromId == next_chrom_id {
            //    next_chrom_id += 1;
            //}
            let (section, _) = BigWig::write_section(file, items_in_section, chromId, current_offset)?;
            sections.push(section);

            for mut zoom in &mut zooms {
                match &zoom.2 {
                    None => (),
                    Some(zoom2) => {
                        assert!(lastchrom == zoom2.chrom);
                        write_zoom(&mut zoom.1, &zoom.2, chromId)?;
                        zoom.3 += 1;
                        zoom.2 = None;
                        println!("File size ({:?}): {:?}", zoom.0, zoom.1.seek(SeekFrom::End(0))?)
                    }
                }
            }
        }

        println!("Zooms: {:?}", zooms);
        let zooms_out: Vec<TempZoomInfo> = zooms.into_iter().map(|z| (z.0, z.1.into_inner().expect("Can't get inner file"), z.3)).collect();
        let final_summary = match summary {
            None => Summary {
                bases_covered: 0,
                min_val: 0.0,
                max_val: 0.0,
                sum: 0.0,
                sum_squares: 0.0,
            },
            Some(summary) => summary,
        };
        Ok((sections, final_summary, zooms_out))
    }

    fn write_zooms(&mut self, zooms: Vec<TempZoomInfo>, zoom_entries: &mut Vec<ZoomHeader>) -> std::io::Result<()> {
        // TODO: param
        const DO_COMPRESS: bool = true;
        const ITEMS_PER_SLOT: u32 = 1024;

        for mut zoom in zooms {
            let (sections, zoom_data_offset, zoom_index_offset) = {
                let mut file = self.get_buf_writer();
                let zoom_data_offset = file.seek(SeekFrom::Current(0))?;
                zoom.1.seek(SeekFrom::Start(0))?;
                let mut sections: Vec<Section> = vec![];
                let mut buf: Vec<u8> = vec![];
                let mut state: Option<(u32, u32, u32)> = None;
                let mut items_in_section: u32 = 0;

                let mut write_section = |buf: Vec<u8>, state: &Option<(u32, u32, u32)>| -> std::io::Result<Section> {
                    let current_offset = file.seek(SeekFrom::Current(0))?;
                    let out_bytes = if DO_COMPRESS {
                        let mut e = ZlibEncoder::new(Vec::new(), Compression::default());
                        e.write_all(&buf)?;
                        e.finish()?
                    } else {
                        buf
                    };
                    file.write_all(&out_bytes)?;
                    let end_offset = file.seek(SeekFrom::Current(0))?;
                    Ok(Section {
                        offset: current_offset,
                        size: end_offset - current_offset,
                        chrom: state.unwrap().0,
                        start: state.unwrap().1,
                        end: state.unwrap().2,
                    })
                };
                for _ in 0..zoom.2 {
                    let chrom = read_u32(&mut zoom.1, false)?;
                    let start = read_u32(&mut zoom.1, false)?;
                    let end = read_u32(&mut zoom.1, false)?;
                    let mut rest = vec![0u8; 20];
                    zoom.1.read_exact(&mut rest)?;
                    //println!("Read zoom entry: {:?} {:?} {:?} {:?}", chrom, start, end, rest);
                    if items_in_section >= ITEMS_PER_SLOT || (state.is_some() && state.as_ref().unwrap().0 != chrom) {
                        let section = write_section(buf, &state)?;
                        sections.push(section);
                        buf = vec![];
                        state = None;
                        items_in_section = 0;
                    }
                    buf.write_u32::<NativeEndian>(chrom)?;
                    buf.write_u32::<NativeEndian>(start)?;
                    buf.write_u32::<NativeEndian>(end)?;
                    buf.write_all(&rest)?;
                    match state {
                        None => {
                            state = Some((chrom, start, end));
                        },
                        Some((_chrom, state_start, _end)) => state = Some((chrom, state_start, end)),
                    }
                    items_in_section += 1;
                }
                let section = write_section(buf, &state)?;
                sections.push(section);
                
                let zoom_index_offset = file.seek(SeekFrom::Current(0))?;
                (sections, zoom_data_offset, zoom_index_offset)
            };
            self.write_rtreeindex(sections)?;

            zoom_entries.push(ZoomHeader {
                reduction_level: zoom.0,
                data_offset: zoom_data_offset,
                index_offset: zoom_index_offset,
            });
        }

        Ok(())
    }

    fn write_rtreeindex(&mut self, sections: Vec<Section>) -> std::io::Result<()> {
        const BLOCK_SIZE: u32 = 256;

        let section_count = sections.len() as u64;

        let mut current_nodes: Vec<RTreeNode> = sections.into_iter().map(|s| RTreeNode {
            start_chrom_idx: s.chrom,
            start_base: s.start,
            end_chrom_idx: s.chrom,
            end_base: s.end,
            kind: RTreeNodeType::Leaf {
                offset: s.offset,
                size: s.size,
            },
        }).collect();
        let mut levels = 0;
        let nodes: RTreeNodeList<RTreeNode> = loop {
            println!("Current_nodes (at level {}): {:?}", levels, current_nodes.len());
            levels += 1;

            if current_nodes.len() < BLOCK_SIZE as usize {
                break RTreeNodeList::<RTreeNode> {
                    nodes: current_nodes
                }
            }

            let mut start_chrom_idx = 0;
            let mut start_base = 0;
            let mut end_chrom_idx = 0;
            let mut end_base = 0;
            let mut next_nodes: Vec<RTreeNode> = vec![];
            let mut current_group: Vec<RTreeNode> = vec![];
            let mut node_iter = current_nodes.into_iter();
            loop {
                let next_node = node_iter.next();
                match next_node {
                    None => {
                        println!("Remaining nodes at complete: {}", current_group.len());
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
                        break
                    },
                    Some(node) => {

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
            current_nodes = next_nodes;
        };
        //println!("Nodes ({:?}): {:?}", nodes.nodes.len(), nodes);
        //println!("Levels: {:?}", levels);


        const NODEHEADER_SIZE: u64 = 1 + 1 + 2;
        const NON_LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8;
        const LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8 + 8;

        let mut file = self.get_buf_writer();

        let mut index_offsets: Vec<u64> = vec![0u64; levels as usize];

        fn calculate_offsets(mut index_offsets: &mut Vec<u64>, trees: &RTreeNodeList<RTreeNode>, level: usize, levels: usize) -> std::io::Result<()> {
            let mut isleaf: i8 = -1;
            let mut current_index_size = NODEHEADER_SIZE;
            for tree in trees.nodes.iter() {
                match &tree.kind {
                    RTreeNodeType::Leaf { .. } => {
                        if level != 0 {
                            panic!("Leaf node found at level {}", level)    
                        }
                        if isleaf == 0 {
                            panic!("Mixed node types at level {}", level);
                        }
                        isleaf = 1;
                        current_index_size += LEAFNODE_SIZE;
                    },
                    RTreeNodeType::NonLeaf { children, .. } => {
                        if level == 0 {
                            panic!("Non leaf node found at level 0")
                        }
                        if isleaf == 1 {
                            panic!("Mixed node types at level {}", level);
                        }
                        isleaf = 0;
                        calculate_offsets(&mut index_offsets, &children, level - 1, levels)?;
                        current_index_size += NON_LEAFNODE_SIZE;
                    },
                }
            }
            index_offsets[level] += current_index_size;
            Ok(())
        }

        if levels > 0 {
            calculate_offsets(&mut index_offsets, &nodes, levels - 1, levels)?;
        }
        println!("index Offsets: {:?}", index_offsets);


        const NON_LEAFNODE_FULL_BLOCK_SIZE: u64 = NODEHEADER_SIZE + NON_LEAFNODE_SIZE * BLOCK_SIZE as u64;
        const LEAFNODE_FULL_BLOCK_SIZE: u64 = NODEHEADER_SIZE + LEAFNODE_SIZE * BLOCK_SIZE as u64;
        fn write_tree(mut file: &mut Write, index_offsets: &Vec<u64>, trees: &RTreeNodeList<RTreeNode>, curr_level: usize, dest_level: usize, mut offset: &mut u64) -> std::io::Result<()> {
            assert!(curr_level >= dest_level);
            if curr_level != dest_level {
                for tree in trees.nodes.iter() {
                    match &tree.kind {
                        RTreeNodeType::Leaf { .. } => panic!("Leaf node found at level {}", curr_level),
                        RTreeNodeType::NonLeaf { children, .. } => {
                            write_tree(&mut file, &index_offsets, &children, curr_level - 1, dest_level, &mut offset)?;
                        },
                    }
                }
                return Ok(())
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
            *offset += index_offsets[curr_level];
            for (idx, node) in trees.nodes.iter().enumerate() {
                file.write_u32::<NativeEndian>(node.start_chrom_idx)?;
                file.write_u32::<NativeEndian>(node.start_base)?;
                file.write_u32::<NativeEndian>(node.end_chrom_idx)?;
                file.write_u32::<NativeEndian>(node.end_base)?;
                match &node.kind {
                    RTreeNodeType::Leaf { offset, size } => {
                        file.write_u64::<NativeEndian>(*offset)?;
                        file.write_u64::<NativeEndian>(*size)?;
                    },
                    RTreeNodeType::NonLeaf { .. } => {
                        let full_size = if curr_level >= 2 {
                            NON_LEAFNODE_FULL_BLOCK_SIZE
                        } else {
                            LEAFNODE_FULL_BLOCK_SIZE
                        };
                        let child_offset: u64 = *offset + idx as u64 * full_size;
                        //println!("Child node offset: {}; Added: {}", child_offset, idx as u64 * full_size);
                        file.write_u64::<NativeEndian>(child_offset)?;
                    },
                }
            }
            Ok(())
        }

        // TODO remove
        const ITEMS_PER_SLOT: u32 = 1024;

        let end_of_data = file.seek(SeekFrom::Current(0))?;
        // TODO: handle case with no data
        {
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
        println!("Start of index: {}", current_offset);
        for level in (0..levels).rev() {
            write_tree(&mut file, &index_offsets, &nodes, levels - 1, level, &mut current_offset)?;
            println!("End of index level {}: {}", level, file.seek(SeekFrom::Current(0))?);
        }

        //assert!(file.seek(SeekFrom::Current(0))? == next_offset);

        Ok(())
    }
}


// TODO: switch to byteordered, which allows for runtime endianness
fn read_u8(f: &mut Read, bigendian: bool) -> std::io::Result<u8> {
    if bigendian {
        return f.read_u8()
    } else {
        return f.read_u8()
    }
}

fn read_u16(f: &mut dyn Read, bigendian: bool) -> std::io::Result<u16> {
    if bigendian {
        return f.read_u16::<BigEndian>()
    } else {
        return f.read_u16::<LittleEndian>()
    }
}

fn read_u32(f: &mut dyn Read, bigendian: bool) -> std::io::Result<u32> {
    if bigendian {
        return f.read_u32::<BigEndian>()
    } else {
        return f.read_u32::<LittleEndian>()
    }
}

fn read_u64(f: &mut dyn Read, bigendian: bool) -> std::io::Result<u64> {
    if bigendian {
        return f.read_u64::<BigEndian>()
    } else {
        return f.read_u64::<LittleEndian>()
    }
}

fn read_f32(f: &mut dyn Read, bigendian: bool) -> std::io::Result<f32> {
    if bigendian {
        return f.read_f32::<BigEndian>()
    } else {
        return f.read_f32::<LittleEndian>()
    }
}