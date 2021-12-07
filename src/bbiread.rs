use std::fs::File;
use std::io::{self, BufReader, Cursor, Read, Seek, SeekFrom};
use std::marker::PhantomData;
use std::vec::Vec;

use byteordered::{ByteOrdered, Endianness};

use crate::bigwig::{
    BBIFile, Summary, ZoomHeader, ZoomRecord, BIGBED_MAGIC, BIGWIG_MAGIC, CHROM_TREE_MAGIC,
    CIR_TREE_MAGIC,
};
use crate::mem_cached_file::MemCachedRead;
use crate::seekableread::SeekableRead;

#[derive(Debug)]
pub struct Block {
    pub offset: u64,
    pub size: u64,
}

#[derive(Clone, Debug)]
pub struct BBIHeader {
    pub endianness: Endianness,

    pub(crate) _version: u16,
    pub(crate) zoom_levels: u16,
    pub(crate) chromosome_tree_offset: u64,
    pub(crate) full_data_offset: u64,
    pub(crate) full_index_offset: u64,
    pub(crate) _field_count: u16,
    pub(crate) _defined_field_count: u16,
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

pub trait BBIRead<R: SeekableRead> {
    fn get_info(&self) -> &BBIFileInfo;

    fn autosql(&mut self) -> io::Result<String>;

    fn ensure_reader(&mut self) -> io::Result<&mut ByteOrdered<BufReader<R>, Endianness>>;

    fn ensure_mem_cached_reader(
        &mut self,
    ) -> io::Result<
        ByteOrdered<BufReader<MemCachedRead<ByteOrdered<BufReader<R>, Endianness>>>, Endianness>,
    >;

    /// Manually close the open file descriptor (if it exists). If any operations are performed after this is called, the file descriptor will be reopened.
    fn close(&mut self);

    fn get_chroms(&self) -> Vec<ChromAndSize>;

    /// This assumes the file is at the cir tree start
    fn search_cir_tree(
        &mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> io::Result<Vec<Block>> {
        // TODO: Move anything relying on self out to separate method
        let chrom_ix = {
            let chrom_info = &self.get_info().chrom_info;
            let chrom = chrom_info.iter().find(|&x| x.name == chrom_name);
            //println!("Chrom: {:?}", chrom);
            match chrom {
                Some(c) => c.id,
                None => {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("{} not found.", chrom_name),
                    ))
                }
            }
        };

        let mut file = self.ensure_mem_cached_reader()?;
        let magic = file.read_u32()?;
        if magic != CIR_TREE_MAGIC {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "Invalid file format: CIR_TREE_MAGIC does not match.",
            ));
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

    fn get_overlapping_blocks(
        &mut self,
        chrom_name: &str,
        start: u32,
        end: u32,
    ) -> io::Result<Vec<Block>> {
        let full_index_offset = self.get_info().header.full_index_offset;

        let file = self.ensure_reader()?;
        file.seek(SeekFrom::Start(full_index_offset))?;

        self.search_cir_tree(chrom_name, start, end)
    }
}

pub(crate) fn read_info<R: SeekableRead>(
    file: BufReader<R>,
) -> Result<BBIFileInfo, BBIFileReadInfoError> {
    let mut file = ByteOrdered::runtime(file, Endianness::Little);

    let magic = file.read_u32()?;
    // println!("Magic {:x?}", magic);
    let filetype = match magic {
        _ if magic == BIGWIG_MAGIC.to_be() => {
            file = file.into_opposite();
            BBIFile::BigWig
        }
        _ if magic == BIGWIG_MAGIC.to_le() => BBIFile::BigWig,
        _ if magic == BIGBED_MAGIC.to_be() => {
            file = file.into_opposite();
            BBIFile::BigBed
        }
        _ if magic == BIGBED_MAGIC.to_le() => BBIFile::BigBed,
        _ => return Err(BBIFileReadInfoError::UnknownMagic),
    };

    let _version = file.read_u16()?;

    // TODO: should probably handle versions < 3
    let zoom_levels = file.read_u16()?;
    let chromosome_tree_offset = file.read_u64()?;
    let full_data_offset = file.read_u64()?;
    let full_index_offset = file.read_u64()?;
    let _field_count = file.read_u16()?;
    let _defined_field_count = file.read_u16()?;
    let auto_sql_offset = file.read_u64()?;
    let total_summary_offset = file.read_u64()?;
    let uncompress_buf_size = file.read_u32()?;
    let _reserved = file.read_u64()?;

    let header = BBIHeader {
        endianness: file.endianness(),
        _version,
        zoom_levels,
        chromosome_tree_offset,
        full_data_offset,
        full_index_offset,
        _field_count,
        _defined_field_count,
        auto_sql_offset,
        total_summary_offset,
        uncompress_buf_size,
    };

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
    assert_eq!(val_size, 8u32);

    let mut chrom_info = Vec::with_capacity(item_count as usize);
    read_chrom_tree_block(&mut file, &mut chrom_info, key_size)?;
    chrom_info.sort_by(|c1, c2| c1.name.cmp(&c2.name));

    let info = BBIFileInfo {
        filetype,
        header,
        zoom_headers,
        chrom_info,
    };

    Ok(info)
}

fn read_zoom_headers<R: SeekableRead>(
    file: &mut ByteOrdered<BufReader<R>, Endianness>,
    header: &BBIHeader,
) -> io::Result<Vec<ZoomHeader>> {
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

fn read_chrom_tree_block<R: SeekableRead>(
    f: &mut ByteOrdered<BufReader<R>, Endianness>,
    chroms: &mut Vec<ChromInfo>,
    key_size: u32,
) -> io::Result<()> {
    let isleaf = f.read_u8()?;
    let _reserved = f.read_u8()?;
    let count = f.read_u16()?;

    if isleaf == 1 {
        for _ in 0..count {
            let mut key_bytes = vec![0u8; key_size as usize];
            f.read_exact(&mut key_bytes)?;
            let key_string = match String::from_utf8(key_bytes) {
                Ok(s) => s.trim_matches(char::from(0)).to_owned(),
                Err(_) => {
                    return Err(io::Error::new(
                        io::ErrorKind::Other,
                        "Invalid file format: Invalid utf-8 string.",
                    ))
                }
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
        // First, go through and get child blocks
        let mut children: Vec<u64> = vec![];
        children.reserve_exact(count as usize);
        for _ in 0..count {
            // We don't need this, but have to read it
            let mut key_bytes = vec![0u8; key_size as usize];
            f.read_exact(&mut key_bytes)?;

            // TODO: could add specific find here by comparing key string
            let child_offset = f.read_u64()?;
            children.push(child_offset);
        }
        // Then go through each child block
        for child in children {
            f.seek(SeekFrom::Start(child))?;
            read_chrom_tree_block(f, chroms, key_size)?;
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
fn overlaps(
    chromq: u32,
    chromq_start: u32,
    chromq_end: u32,
    chromb1: u32,
    chromb1_start: u32,
    chromb2: u32,
    chromb2_end: u32,
) -> bool {
    compare_position(chromq, chromq_start, chromb2, chromb2_end) <= 0
        && compare_position(chromq, chromq_end, chromb1, chromb1_start) >= 0
}

pub(crate) fn search_overlapping_blocks<R: SeekableRead>(
    file: &mut ByteOrdered<R, Endianness>,
    chrom_ix: u32,
    start: u32,
    end: u32,
    blocks: &mut Vec<Block>,
) -> io::Result<()> {
    //println!("Searching for overlapping blocks at {:?}. Searching {:?}:{:?}-{:?}", self.current_file_offset()?, chrom_ix, start, end);

    let isleaf: u8 = file.read_u8()?;
    assert!(isleaf == 1 || isleaf == 0, "Unexpected isleaf: {}", isleaf);
    let _reserved = file.read_u8()?;
    let count: u16 = file.read_u16()?;
    //eprintln!("Index: {:?} {:?} {:?}", isleaf, _reserved, count);

    let mut childblocks: Vec<u64> = vec![];
    for _ in 0..count {
        let start_chrom_ix = file.read_u32()?;
        let start_base = file.read_u32()?;
        let end_chrom_ix = file.read_u32()?;
        let end_base = file.read_u32()?;
        let block_overlaps = overlaps(
            chrom_ix,
            start,
            end,
            start_chrom_ix,
            start_base,
            end_chrom_ix,
            end_base,
        );
        if isleaf == 1 {
            let data_offset = file.read_u64()?;
            let data_size = file.read_u64()?;
            if block_overlaps {
                //eprintln!("Overlaps (leaf): {:?}:{:?}-{:?} with {:?}:{:?}-{:?}:{:?} {:?} {:?}", chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset, data_size);
                blocks.push(Block {
                    offset: data_offset,
                    size: data_size,
                });
            }
        } else {
            let data_offset = file.read_u64()?;
            if block_overlaps {
                //eprintln!("Overlaps (non-leaf): {:?}:{:?}-{:?} with {:?}:{:?}-{:?}:{:?} {:?}", chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base, data_offset);
                childblocks.push(data_offset);
            }
        }
    }
    for childblock in childblocks {
        //eprintln!("Seeking to {:?}", childblock);
        file.seek(SeekFrom::Start(childblock))?;
        search_overlapping_blocks(file, chrom_ix, start, end, blocks)?;
    }
    Ok(())
}

pub fn get_filetype(path: &str) -> io::Result<Option<BBIFile>> {
    let mut file = ByteOrdered::runtime(File::open(path)?, Endianness::Little);
    let magic = file.read_u32()?;
    let file_type = match magic {
        _ if magic == BIGWIG_MAGIC.to_be() => Some(BBIFile::BigWig),
        _ if magic == BIGWIG_MAGIC.to_le() => Some(BBIFile::BigWig),
        _ if magic == BIGBED_MAGIC.to_be() => Some(BBIFile::BigBed),
        _ if magic == BIGBED_MAGIC.to_le() => Some(BBIFile::BigBed),
        _ => None,
    };
    Ok(file_type)
}

/// Gets the data (uncompressed, if applicable) from a given block
pub(crate) fn get_block_data<S: SeekableRead, B: BBIRead<S>>(
    bbifile: &mut B,
    block: &Block,
    known_offset: u64,
) -> io::Result<ByteOrdered<Cursor<Vec<u8>>, Endianness>> {
    use libdeflater::Decompressor;

    let (endianness, uncompress_buf_size) = {
        let info = bbifile.get_info();
        (
            info.header.endianness,
            info.header.uncompress_buf_size as usize,
        )
    };
    let file = bbifile.ensure_reader()?;

    // TODO: Could minimize this by chunking block reads
    if known_offset != block.offset {
        file.seek(SeekFrom::Start(block.offset))?;
    }

    let mut raw_data = vec![0u8; block.size as usize];
    file.read_exact(&mut raw_data)?;
    let block_data: Vec<u8> = if uncompress_buf_size > 0 {
        let mut decompressor = Decompressor::new();
        let mut outbuf = vec![0; uncompress_buf_size];
        let decompressed = decompressor
            .zlib_decompress(&raw_data, &mut outbuf)
            .unwrap();
        outbuf.truncate(decompressed);
        outbuf
    } else {
        raw_data
    };

    Ok(ByteOrdered::runtime(Cursor::new(block_data), endianness))
}

pub(crate) fn get_zoom_block_values<S: SeekableRead, B: BBIRead<S>>(
    bbifile: &mut B,
    block: Block,
    known_offset: &mut u64,
    start: u32,
    end: u32,
) -> io::Result<Box<dyn Iterator<Item = ZoomRecord> + Send>> {
    let mut data_mut = get_block_data(bbifile, &block, *known_offset)?;
    let len = data_mut.inner_mut().get_mut().len();
    assert_eq!(len % (4 * 8), 0);
    let itemcount = len / (4 * 8);
    let mut records = Vec::with_capacity(itemcount);

    for _ in 0..itemcount {
        let chrom_id = data_mut.read_u32()?;
        let chrom_start = data_mut.read_u32()?;
        let chrom_end = data_mut.read_u32()?;
        let bases_covered = u64::from(data_mut.read_u32()?);
        let min_val = f64::from(data_mut.read_f32()?);
        let max_val = f64::from(data_mut.read_f32()?);
        let sum = f64::from(data_mut.read_f32()?);
        let sum_squares = f64::from(data_mut.read_f32()?);
        if chrom_end >= start && chrom_start <= end {
            records.push(ZoomRecord {
                chrom: chrom_id,
                start: chrom_start,
                end: chrom_end,
                summary: Summary {
                    total_items: 0,
                    bases_covered,
                    min_val,
                    max_val,
                    sum,
                    sum_squares,
                },
            });
        }
    }

    *known_offset = block.offset + block.size;
    Ok(Box::new(records.into_iter()))
}

pub(crate) struct ZoomIntervalIter<'a, I, S, B>
where
    I: Iterator<Item = Block> + Send,
    S: SeekableRead,
    B: BBIRead<S>,
{
    bbifile: &'a mut B,
    known_offset: u64,
    blocks: I,
    vals: Option<Box<dyn Iterator<Item = ZoomRecord> + Send + 'a>>,
    start: u32,
    end: u32,
    _phantom: PhantomData<S>,
}

impl<'a, I, S, B> ZoomIntervalIter<'a, I, S, B>
where
    I: Iterator<Item = Block> + Send,
    S: SeekableRead,
    B: BBIRead<S>,
{
    pub fn new(bbifile: &'a mut B, blocks: I, start: u32, end: u32) -> Self {
        ZoomIntervalIter {
            bbifile,
            known_offset: 0,
            blocks,
            vals: None,
            start,
            end,
            _phantom: PhantomData,
        }
    }
}

impl<'a, I, S, B> Iterator for ZoomIntervalIter<'a, I, S, B>
where
    I: Iterator<Item = Block> + Send,
    S: SeekableRead,
    B: BBIRead<S>,
{
    type Item = io::Result<ZoomRecord>;

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
                    let current_block = self.blocks.next()?;
                    match get_zoom_block_values(
                        self.bbifile,
                        current_block,
                        &mut self.known_offset,
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
