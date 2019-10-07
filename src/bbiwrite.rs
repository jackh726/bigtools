use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::pin::Pin;

use byteorder::{NativeEndian, WriteBytesExt};

use flate2::Compression;
use flate2::write::ZlibEncoder;

use futures::try_join;
use futures::channel::mpsc::{channel, Receiver};
use futures::executor::{ThreadPool};
use futures::future::{Either, Future, FutureExt};
use futures::stream::StreamExt;
use futures::task::SpawnExt;

use serde::{Serialize, Deserialize};

use crate::filebufferedchannel;
use crate::idmap::IdMap;
use crate::tell::Tell;
use crate::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};

use crate::bigwig::{CHROM_TREE_MAGIC, CIR_TREE_MAGIC, ZoomHeader, Summary, ZoomRecord};


pub(crate) struct ZoomInfo {
    pub(crate) resolution: u32,
    pub(crate) data: File,
    pub(crate) sections: Box<dyn Iterator<Item=Section>>,
}

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
    pub zoom_sizes: Vec<u32>,
}

pub(crate) const DEFAULT_ZOOM_SIZES: [u32; 11] = [10, 40, 160, 640, 2_560, 10_240, 40_960, 163_840, 655_360, 2_621_440, 10_485_760];

#[derive(Debug)]
pub enum WriteGroupsError {
    InvalidInput(String),
    IoError(io::Error),
}

impl std::fmt::Display for WriteGroupsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(message) => write!(f, "{}", message),
            Self::IoError(e) => write!(f, "{}", e),
        }
    }
}

impl From<io::Error> for WriteGroupsError {
    fn from(error: io::Error) -> Self {
        WriteGroupsError::IoError(error)
    }
}

impl<W> From<io::IntoInnerError<W>> for WriteGroupsError {
    fn from(error: io::IntoInnerError<W>) -> Self {
        WriteGroupsError::IoError(error.into())
    }
}

pub struct TempZoomInfo {
    pub resolution: u32,
    pub data_write_future: Box<dyn Future<Output=Result<usize, WriteGroupsError>> + Send + Unpin>,
    pub data: TempFileBuffer<TempFileBufferWriter<File>>,
    pub sections: filebufferedchannel::Receiver<Section>,   
}

pub(crate) struct ChromProcessingInput {
    pub(crate) zooms_channels: Vec<futures::channel::mpsc::Sender<Pin<Box<dyn Future<Output=io::Result<SectionData>> + Send>>>>,
    pub(crate) ftx: futures::channel::mpsc::Sender<Pin<Box<dyn Future<Output=io::Result<SectionData>> + Send>>>,
}

pub struct ChromProcessingOutput {
    pub sections: filebufferedchannel::Receiver<Section>,
    pub data: TempFileBuffer<File>,
    pub data_write_future: Box<dyn Future<Output=Result<usize, WriteGroupsError>> + Send + Unpin>,
    pub zooms: Vec<TempZoomInfo>,
}

pub struct ChromGroupRead {
    pub summary_future: Pin<Box<dyn Future<Output=Result<Summary, WriteGroupsError>> + Send>>,
    pub processing_output: ChromProcessingOutput,
}

pub trait ChromGroupReadStreamingIterator {
    fn next(&mut self) -> Result<Option<Either<ChromGroupRead, IdMap<String>>>, WriteGroupsError>;
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

pub(crate) async fn encode_zoom_section(compress: bool, items_in_section: Vec<ZoomRecord>) -> io::Result<SectionData> {
    let mut bytes: Vec<u8> = vec![];

    let start = items_in_section[0].start;
    let end = items_in_section[items_in_section.len() - 1].end;

    let chrom = items_in_section[0].chrom;
    for item in items_in_section.iter() {
        bytes.write_u32::<NativeEndian>(item.chrom)?;
        bytes.write_u32::<NativeEndian>(item.start)?;
        bytes.write_u32::<NativeEndian>(item.end)?;
        bytes.write_u32::<NativeEndian>(item.summary.bases_covered as u32)?;
        bytes.write_f32::<NativeEndian>(item.summary.min_val as f32)?;
        bytes.write_f32::<NativeEndian>(item.summary.max_val as f32)?;
        bytes.write_f32::<NativeEndian>(item.summary.sum as f32)?;
        bytes.write_f32::<NativeEndian>(item.summary.sum_squares as f32)?; 
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
    let non_leafnode_full_block_size: u64 = NODEHEADER_SIZE + NON_LEAFNODE_SIZE * u64::from(options.block_size);
    let leafnode_full_block_size: u64 = NODEHEADER_SIZE + LEAFNODE_SIZE * u64::from(options.block_size);
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
        let mut zoom_file = zoom.data;
        let zoom_size = zoom_file.seek(SeekFrom::End(0))?;
        if zoom_size > (data_size / 2) {
            continue;
        }
        let zoom_data_offset = file.tell()?;

        let mut current_offset = zoom_data_offset;
        let sections_iter = zoom.sections.map(|mut section| {
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
            reduction_level: zoom.resolution,
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

pub(crate) async fn write_vals<V>(
    mut vals_iter: V,
    file: BufWriter<File>,
    options: &BBIWriteOptions,
)
-> Result<(
    IdMap<String>,
    Summary,
    BufWriter<File>,
    Box<dyn Iterator<Item=Section> + 'static>,
    Vec<ZoomInfo>
), WriteGroupsError>
where V : ChromGroupReadStreamingIterator + Send {
    // Zooms have to be double-buffered: first because chroms could be processed in parallel and second because we don't know the offset of each zoom immediately
    type ZoomValue = (Vec<Box<dyn Iterator<Item=Section>>>, TempFileBuffer<File>, Option<TempFileBufferWriter<File>>);
    let mut zooms_map: HashMap<u32, ZoomValue> = options.zoom_sizes.iter().map(|size| -> io::Result<_> {
        let section_iter: Vec<Box<dyn Iterator<Item=Section>>> = vec![];
        let (buf, write): (TempFileBuffer<File>, TempFileBufferWriter<File>) = TempFileBuffer::new()?;
        let value = (section_iter, buf, Some(write));
        Ok((*size, value))
    }).collect::<io::Result<_>>()?;

    let mut section_iter: Vec<Box<dyn Iterator<Item=Section>>> = vec![];
    let mut raw_file = file.into_inner()?;

    let mut summary: Option<Summary> = None;

    let chrom_ids = loop {
        let next = vals_iter.next()?;
        match next {
            Some(Either::Left(read)) => {
                let ChromGroupRead { summary_future, processing_output: ChromProcessingOutput { sections, mut data, data_write_future, mut zooms } } = read;
                // If we concurrently processing multiple chromosomes, the section buffer might have written some or all to a separate file
                // Switch that processing output to the real file
                data.switch(raw_file);
                for TempZoomInfo { resolution: size, data: buf, .. } in zooms.iter_mut() {
                    let zoom = zooms_map.get_mut(size).unwrap();
                    let writer = zoom.2.take().unwrap();
                    buf.switch(writer);
                }
                
                // All the futures are actually just handles, so these are purely for the result
                let (chrom_summary, _num_sections) = try_join!(summary_future, data_write_future)?;
                section_iter.push(Box::new(sections.into_iter()));
                raw_file = data.await_real_file();

                for TempZoomInfo { resolution, data_write_future, data, sections } in zooms.into_iter() {
                    let zoom = zooms_map.get_mut(&resolution).unwrap();
                    let _num_sections = data_write_future.await?;
                    zoom.0.push(Box::new(sections.into_iter()));
                    zoom.2.replace(data.await_real_file());
                }

                match &mut summary {
                    None => {
                        summary = Some(chrom_summary)
                    },
                    Some(summary) => {
                        summary.total_items += chrom_summary.total_items;
                        summary.bases_covered += chrom_summary.bases_covered;
                        summary.min_val = summary.min_val.min(chrom_summary.min_val);
                        summary.max_val = summary.max_val.max(chrom_summary.max_val);
                        summary.sum += chrom_summary.sum;
                        summary.sum_squares += chrom_summary.sum_squares;
                    }
                }
            },
            Some(Either::Right(chrom_ids)) => break chrom_ids,
            None => unreachable!(),
        }
    };

    let summary_complete = summary.unwrap_or(Summary {
        total_items: 0,
        bases_covered: 0,
        min_val: 0.0,
        max_val: 0.0,
        sum: 0.0,
        sum_squares: 0.0,
    });

    let zoom_infos: Vec<ZoomInfo> = zooms_map.into_iter().map(|(size, zoom)| {
        drop(zoom.2);
        let zoom_iter: Box<dyn Iterator<Item=Section> + 'static> = Box::new(zoom.0.into_iter().flat_map(|s| s));
        let closed_file = zoom.1.await_temp_file();
        ZoomInfo {
            resolution: size,
            data: closed_file,
            sections: zoom_iter,
        }
    }).collect();
    let section_iter: Box<dyn Iterator<Item=Section> + 'static> = Box::new(section_iter.into_iter().flat_map(|s| s));
    Ok((chrom_ids, summary_complete, BufWriter::new(raw_file), section_iter, zoom_infos))
}

pub(crate) async fn write_data<W: Write>(
    mut data_file: W,
    mut section_sender: filebufferedchannel::Sender<Section>,
    mut frx: Receiver<impl Future<Output=io::Result<SectionData>> + Send >
) -> Result<usize, WriteGroupsError> {
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
            size,
        }).expect("Couldn't send section.");
        current_offset += size;
    }
    Ok(total)
}

pub(crate) fn get_chromprocessing(
    pool: &mut ThreadPool,
    options: &BBIWriteOptions,
) -> io::Result<(ChromProcessingInput, ChromProcessingOutput)> {
    let (ftx, frx) = channel::<_>(100);

    let (sections_handle, buf, section_receiver) = {
        let (buf, write) = TempFileBuffer::new()?;
        let file = BufWriter::new(write);

        let (section_sender, section_receiver) = filebufferedchannel::channel::<Section>(200);
        let (sections_remote, sections_handle) = write_data(file, section_sender, frx).remote_handle();
        pool.spawn(sections_remote).expect("Couldn't spawn future.");
        (sections_handle, buf, section_receiver)
    };

    let processed_zooms: Vec<_> = options.zoom_sizes.iter().map(|size| -> io::Result<_> {
        let (ftx, frx) = channel::<_>(100);
        let (buf, write) = TempFileBuffer::new()?;
        let file = BufWriter::new(write);

        let (section_sender, section_receiver) = filebufferedchannel::channel::<Section>(200);
        let (remote, handle) = write_data(file, section_sender, frx).remote_handle();
        pool.spawn(remote).expect("Couldn't spawn future.");
        let zoom_info = TempZoomInfo {
            resolution: *size,
            data_write_future: Box::new(handle),
            data: buf,
            sections: section_receiver,
        };
        Ok((zoom_info, ftx))
    }).collect::<io::Result<_>>()?;
    let (zoom_infos, zooms_channels): (Vec<_>, Vec<_>) = processed_zooms.into_iter().unzip();

    Ok((
        ChromProcessingInput {
            zooms_channels,
            ftx,
        },
        ChromProcessingOutput {
            sections: section_receiver,
            data: buf,
            data_write_future: Box::new(sections_handle),
            zooms: zoom_infos,
        }
    ))
}