use std::io::{self, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::fs::File;

use byteorder::{NativeEndian, WriteBytesExt};

use flate2::Compression;
use flate2::write::ZlibEncoder;

use serde::{Serialize, Deserialize};

use crate::tell::Tell;

use crate::bigwig::{CHROM_TREE_MAGIC, CIR_TREE_MAGIC, Value, ZoomHeader};


pub(crate) type ZoomInfo = (u32 /* resolution */, File /* Temp file that contains data */, Box<dyn Iterator<Item=Section>> /* sections */);

#[derive(Debug)]
pub(crate) struct SectionData {
    pub(crate) chrom: u32,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) data: Vec<u8>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Section {
    pub(crate) chrom: u32,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) offset: u64,
    pub(crate) size: u64,
}

#[derive(Debug)]
pub(crate) struct ZoomRecord {
    pub(crate) chrom: u32,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) valid_count: u32,
    pub(crate) min_value: f32,
    pub(crate) max_value: f32,
    pub(crate) sum: f32,
    pub(crate) sum_squares: f32,
}

#[derive(Debug)]
pub(crate) struct RTreeNode {
    pub(crate) start_chrom_idx: u32,
    pub(crate) start_base: u32,
    pub(crate) end_chrom_idx: u32,
    pub(crate) end_base: u32,
    pub(crate) children: RTreeChildren,
}

#[derive(Debug)]
pub(crate) enum RTreeChildren {
    DataSections(Vec<Section>),
    Nodes(Vec<RTreeNode>),
}

#[derive(Clone)]
pub struct BBIWriteOptions {
    pub compress: bool,
    pub items_per_slot: u32,
    pub block_size: u32,
}

pub trait BBIWrite {
}

pub(crate) const MAX_ZOOM_LEVELS: usize = 10;

pub(crate) fn write_blank_headers(file: &mut BufWriter<File>) -> io::Result<()> {
    file.seek(SeekFrom::Start(0))?;
    // Common header
    file.write_all(&[0; 64])?;
    // Zoom levels
    file.write_all(&[0; MAX_ZOOM_LEVELS * 24])?;

    Ok(())
}

pub(crate) fn write_chrom_tree(file: &mut BufWriter<File>, chrom_sizes: std::collections::HashMap<String, u32>, chrom_ids: &std::collections::HashMap<String, u32>) -> io::Result<()> {
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

pub(crate) async fn encode_section(compress: bool, items_in_section: Vec<Value>, chrom_id: u32) -> io::Result<SectionData> {
    let mut bytes: Vec<u8> = vec![];

    let start = items_in_section[0].start;
    let end = items_in_section[items_in_section.len() - 1].end;
    bytes.write_u32::<NativeEndian>(chrom_id)?;
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
        chrom: chrom_id,
        start,
        end,
        data: out_bytes,
    })
}

pub(crate) async fn encode_zoom_section(compress: bool, items_in_section: Vec<ZoomRecord>) -> io::Result<SectionData> {
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

// TODO: it would be cool to output as an iterator so we don't have to store the index in memory
pub(crate) fn get_rtreeindex<S>(sections_stream: S, options: &BBIWriteOptions) -> (RTreeChildren, usize, u64) where S : Iterator<Item=Section> {
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


const NODEHEADER_SIZE: u64 = 1 + 1 + 2;
const NON_LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8;
const LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8 + 8;

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

fn write_tree(
    mut file: &mut BufWriter<File>,
    nodes: &RTreeChildren,
    curr_level: usize,
    dest_level: usize,
    childnode_offset: u64,
    options: &BBIWriteOptions
) -> io::Result<u64> {
    let non_leafnode_full_block_size: u64 = NODEHEADER_SIZE + NON_LEAFNODE_SIZE * options.block_size as u64;
    let leafnode_full_block_size: u64 = NODEHEADER_SIZE + LEAFNODE_SIZE * options.block_size as u64;
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
                non_leafnode_full_block_size
            } else {
                leafnode_full_block_size
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

pub(crate) fn write_rtreeindex(
    file: &mut BufWriter<File>,
    nodes: RTreeChildren,
    levels: usize,
    section_count: u64,
    options: &BBIWriteOptions
) -> io::Result<()> {
    let mut index_offsets: Vec<u64> = vec![0u64; levels as usize];

    calculate_offsets(&mut index_offsets, &nodes, levels);

    let end_of_data = file.seek(SeekFrom::Current(0))?;
    {
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
    for level in (0..=levels).rev() {
        if level > 0 {
            next_offset += index_offsets[level - 1];
        }
        write_tree(file, &nodes, levels, level, next_offset, options)?;
    }

    Ok(())
}

pub(crate) fn write_zooms(mut file: &mut BufWriter<File>, zooms: Vec<ZoomInfo>, data_size: u64, options: &BBIWriteOptions) -> io::Result<Vec<ZoomHeader>> {
    let mut zoom_entries: Vec<ZoomHeader> = vec![];
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

        let (nodes, levels, total_sections) = get_rtreeindex(sections_iter, options);
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
        write_rtreeindex(&mut file, nodes, levels, total_sections, options)?;

        zoom_entries.push(ZoomHeader {
            reduction_level: zoom.0,
            data_offset: zoom_data_offset,
            index_offset: zoom_index_offset,
        });

        zoom_count += 1;
        if zoom_count >= MAX_ZOOM_LEVELS {
            break;
        }
    }

    Ok(zoom_entries)
}
