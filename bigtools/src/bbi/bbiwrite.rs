use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::iter::Flatten;
use std::pin::Pin;
use std::vec;

use byteorder::{NativeEndian, WriteBytesExt};
use crossbeam_channel::unbounded;
use thiserror::Error;

use futures::channel::mpsc as futures_mpsc;
use futures::channel::mpsc::channel;
use futures::future::{Future, FutureExt};
use futures::stream::StreamExt;

use serde::{Deserialize, Serialize};
use tokio::runtime::{Handle, Runtime};

use crate::utils::chromvalues::ChromValues;
use crate::utils::idmap::IdMap;
use crate::utils::tell::Tell;
use crate::utils::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};

use crate::bbi::{Summary, ZoomHeader, ZoomRecord, CHROM_TREE_MAGIC, CIR_TREE_MAGIC};

pub(crate) struct ZoomInfo {
    resolution: u32,
    data: TempFileBuffer<File>,
    sections: Flatten<vec::IntoIter<crossbeam_channel::IntoIter<Section>>>,
}

#[derive(Debug)]
pub(crate) struct SectionData {
    pub(crate) chrom: u32,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) data: Vec<u8>,
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub(crate) struct Section {
    pub(crate) chrom: u32,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) offset: u64,
    pub(crate) size: u64,
}

#[derive(Debug)]
pub(crate) struct RTreeNode {
    start_chrom_idx: u32,
    start_base: u32,
    end_chrom_idx: u32,
    end_base: u32,
    children: RTreeChildren,
}

#[derive(Debug)]
pub(crate) enum RTreeChildren {
    DataSections(Vec<Section>),
    Nodes(Vec<RTreeNode>),
}

/// Options for required input sort order of values
#[derive(Copy, Clone)]
pub enum InputSortType {
    /// Both the chromosomes and start values must be sorted
    ALL,
    /// Start values within a chromosome must be sorted, but chromosomes may be out of order
    START,
    // TODO
    //NONE,
}

/// The default block size used when writing a bbi file
pub const DEFAULT_BLOCK_SIZE: u32 = 256;
/// The default items per slot used when writing a bbi file
pub const DEFAULT_ITEMS_PER_SLOT: u32 = 1024;

/// Options for writing to a bbi file
#[derive(Copy, Clone)]
pub struct BBIWriteOptions {
    pub compress: bool,
    pub items_per_slot: u32,
    pub block_size: u32,
    pub initial_zoom_size: u32,
    pub max_zooms: u32,
    pub input_sort_type: InputSortType,
    pub channel_size: usize,
    pub inmemory: bool,
}

impl Default for BBIWriteOptions {
    fn default() -> Self {
        BBIWriteOptions {
            compress: true,
            items_per_slot: DEFAULT_ITEMS_PER_SLOT,
            block_size: DEFAULT_BLOCK_SIZE,
            initial_zoom_size: 160,
            max_zooms: 10,
            input_sort_type: InputSortType::ALL,
            channel_size: 100,
            inmemory: false,
        }
    }
}

/// Possible errors encountered when processing a chromosome when writing a bbi file
#[derive(Error, Debug)]
pub enum ProcessChromError<SourceError: Error> {
    #[error("{}", .0)]
    InvalidInput(String),
    #[error("{}", .0)]
    InvalidChromosome(String),
    #[error("{}", .0)]
    IoError(#[from] io::Error),
    #[error("{}", .0)]
    SourceError(SourceError),
}

pub(crate) struct TempZoomInfo<SourceError: Error> {
    pub resolution: u32,
    pub data_write_future: Box<
        dyn Future<Output = Result<(usize, usize), ProcessChromError<SourceError>>> + Send + Unpin,
    >,
    pub data: TempFileBuffer<TempFileBufferWriter<File>>,
    pub sections: crossbeam_channel::Receiver<Section>,
}

pub(crate) type ChromProcessingInputSectionChannel = futures::channel::mpsc::Sender<
    Pin<Box<dyn Future<Output = io::Result<(SectionData, usize)>> + Send>>,
>;

const MAX_ZOOM_LEVELS: usize = 10;

pub(crate) fn write_blank_headers(file: &mut BufWriter<File>) -> io::Result<()> {
    file.seek(SeekFrom::Start(0))?;
    // Common header
    file.write_all(&[0; 64])?;
    // Zoom levels
    file.write_all(&[0; MAX_ZOOM_LEVELS * 24])?;

    Ok(())
}

pub(crate) fn write_info<Err: Error>(
    file: &mut BufWriter<File>,
    magic: u32,
    num_zooms: u16,
    chrom_index_start: u64,
    full_data_offset: u64,
    index_start: u64,
    field_count: u16,
    defined_field_count: u16,
    auto_sql_offset: u64,
    total_summary_offset: u64,
    uncompress_buf_size: usize,
    zoom_entries: Vec<ZoomHeader>,
    summary: Summary,
    data_count: u64,
) -> Result<(), ProcessChromError<Err>> {
    file.seek(SeekFrom::Start(0))?;
    file.write_u32::<NativeEndian>(magic)?;
    file.write_u16::<NativeEndian>(4)?;
    file.write_u16::<NativeEndian>(num_zooms)?;
    file.write_u64::<NativeEndian>(chrom_index_start)?;
    file.write_u64::<NativeEndian>(full_data_offset)?;
    file.write_u64::<NativeEndian>(index_start)?;
    file.write_u16::<NativeEndian>(field_count)?; // fieldCount
    file.write_u16::<NativeEndian>(defined_field_count)?; // definedFieldCount
    file.write_u64::<NativeEndian>(auto_sql_offset)?; // autoSQLOffset
    file.write_u64::<NativeEndian>(total_summary_offset)?;
    file.write_u32::<NativeEndian>(uncompress_buf_size as u32)?;
    file.write_u64::<NativeEndian>(0)?; // reserved

    debug_assert!(file.seek(SeekFrom::Current(0))? == 64);

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

    file.seek(SeekFrom::Start(full_data_offset))?;
    file.write_u64::<NativeEndian>(data_count)?;

    file.seek(SeekFrom::End(0))?;
    file.write_u32::<NativeEndian>(magic)?;

    Ok(())
}

pub(crate) fn write_chrom_tree(
    file: &mut BufWriter<File>,
    chrom_sizes: std::collections::HashMap<String, u32>,
    chrom_ids: &std::collections::HashMap<String, u32>,
) -> io::Result<()> {
    let mut chroms: Vec<(&String, &u32)> = chrom_ids.iter().collect();
    chroms.sort_by_key(|v| *v.1);
    //println!("Used chroms {:?}", chroms);

    let item_count = chroms.len() as u64;
    // TODO: for now, always just use the length of chroms (if less than 256). This means we don't have to implement writing non-leaf nodes for now...
    // TODO: make this configurable
    let block_size = std::cmp::max(256, item_count) as u32;
    let max_bytes = chroms
        .iter()
        .map(|a| a.0.as_bytes().len() as u32)
        .fold(0, u32::max);

    file.write_u32::<NativeEndian>(CHROM_TREE_MAGIC)?;
    file.write_u32::<NativeEndian>(block_size)?;
    file.write_u32::<NativeEndian>(max_bytes)?;
    file.write_u32::<NativeEndian>(8)?; // size of Id (u32) + Size (u32)
    file.write_u64::<NativeEndian>(item_count)?;
    file.write_u64::<NativeEndian>(0)?; // Reserved

    // Assuming this is all one block right now
    // TODO: add non-leaf nodes and split blocks
    file.write_u8(1)?;
    file.write_u8(0)?;
    file.write_u16::<NativeEndian>(item_count as u16)?;
    for (chrom, id) in chroms {
        let key_bytes = &mut vec![0u8; max_bytes as usize];
        let chrom_bytes = chrom.as_bytes();
        key_bytes[..chrom_bytes.len()].copy_from_slice(chrom_bytes);
        file.write_all(key_bytes)?;
        file.write_u32::<NativeEndian>(*id)?;
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

pub(crate) async fn encode_zoom_section(
    compress: bool,
    items_in_section: Vec<ZoomRecord>,
) -> io::Result<(SectionData, usize)> {
    use libdeflater::{CompressionLvl, Compressor};

    let mut bytes = Vec::with_capacity(items_in_section.len() * 32);

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

    let (out_bytes, uncompressed_buf_size) = if compress {
        let mut compressor = Compressor::new(CompressionLvl::default());
        let max_sz = compressor.zlib_compress_bound(bytes.len());
        let mut compressed_data = vec![0; max_sz];
        let actual_sz = compressor
            .zlib_compress(&bytes, &mut compressed_data)
            .unwrap();
        compressed_data.resize(actual_sz, 0);
        (compressed_data, bytes.len())
    } else {
        (bytes, 0)
    };

    Ok((
        SectionData {
            chrom,
            start,
            end,
            data: out_bytes,
        },
        uncompressed_buf_size,
    ))
}

// TODO: it would be cool to output as an iterator so we don't have to store the index in memory
pub(crate) fn get_rtreeindex<S>(
    sections_stream: S,
    options: BBIWriteOptions,
) -> (RTreeChildren, usize, u64)
where
    S: Iterator<Item = Section>,
{
    use itertools::Itertools;

    let block_size = options.block_size as usize;
    let mut total_sections = 0;

    let chunks = sections_stream
        .inspect(|_| total_sections += 1)
        .chunks(block_size);
    let mut current_nodes: Vec<RTreeChildren> = chunks
        .into_iter()
        .map(|chunk| {
            let current_chunk: Vec<_> = chunk.collect();
            RTreeChildren::DataSections(current_chunk)
        })
        .collect();
    let mut levels = 0;
    let nodes: RTreeChildren = loop {
        if current_nodes.len() == 1 {
            break current_nodes.pop().unwrap();
        }
        levels += 1;
        let chunks = current_nodes.into_iter().chunks(block_size);
        current_nodes = chunks
            .into_iter()
            .map(|chunk| {
                RTreeChildren::Nodes(
                    chunk
                        .map(|c| match &c {
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
                        })
                        .collect(),
                )
            })
            .collect()
    };

    //println!("Total sections: {:?}", total_sections);
    //println!("Nodes ({:?}): {:?}", nodes.nodes.len(), nodes);
    //println!("Levels: {:?}", levels);
    (nodes, levels, total_sections)
}

const NODEHEADER_SIZE: u64 = 1 + 1 + 2;
const NON_LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8;
const LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8 + 8;

fn calculate_offsets(index_offsets: &mut Vec<u64>, nodes: &RTreeChildren, level: usize) {
    match nodes {
        RTreeChildren::DataSections(_) => (),
        RTreeChildren::Nodes(children) => {
            index_offsets[level - 1] += NODEHEADER_SIZE;
            for child in children {
                index_offsets[level - 1] += NON_LEAFNODE_SIZE;
                calculate_offsets(index_offsets, &child.children, level - 1);
            }
        }
    }
}

fn write_tree<W: Write>(
    file: &mut W,
    nodes: &RTreeChildren,
    curr_level: usize,
    dest_level: usize,
    childnode_offset: u64,
    options: BBIWriteOptions,
) -> io::Result<u64> {
    let non_leafnode_full_block_size: u64 =
        NODEHEADER_SIZE + NON_LEAFNODE_SIZE * u64::from(options.block_size);
    let leafnode_full_block_size: u64 =
        NODEHEADER_SIZE + LEAFNODE_SIZE * u64::from(options.block_size);
    debug_assert!(curr_level >= dest_level);
    if curr_level != dest_level {
        let mut next_offset_offset = 0;
        match nodes {
            RTreeChildren::DataSections(_) => {
                panic!("Datasections found at level: {:?}", curr_level)
            }
            RTreeChildren::Nodes(children) => {
                for child in children {
                    let size = write_tree(
                        file,
                        &child.children,
                        curr_level - 1,
                        dest_level,
                        childnode_offset + next_offset_offset,
                        options,
                    )?;
                    next_offset_offset += size;
                }
            }
        }
        return Ok(next_offset_offset);
    }

    match &nodes {
        RTreeChildren::DataSections(sections) => {
            file.write_u8(1)?;
            file.write_u8(0)?;
            file.write_u16::<NativeEndian>(sections.len() as u16)?;
            for section in sections {
                file.write_u32::<NativeEndian>(section.chrom)?;
                file.write_u32::<NativeEndian>(section.start)?;
                file.write_u32::<NativeEndian>(section.chrom)?;
                file.write_u32::<NativeEndian>(section.end)?;
                file.write_u64::<NativeEndian>(section.offset)?;
                file.write_u64::<NativeEndian>(section.size)?;
            }
            Ok(4 + sections.len() as u64 * 32)
        }
        RTreeChildren::Nodes(children) => {
            file.write_u8(0)?;
            file.write_u8(0)?;
            file.write_u16::<NativeEndian>(children.len() as u16)?;
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
            }
            Ok(children.len() as u64 * full_size)
        }
    }
}

pub(crate) fn write_rtreeindex<W: Write + Seek>(
    file: &mut W,
    nodes: RTreeChildren,
    levels: usize,
    section_count: u64,
    options: BBIWriteOptions,
) -> io::Result<()> {
    let mut index_offsets: Vec<u64> = vec![0u64; levels as usize];

    calculate_offsets(&mut index_offsets, &nodes, levels);

    let end_of_data = file.tell()?;
    file.write_u32::<NativeEndian>(CIR_TREE_MAGIC)?;
    file.write_u32::<NativeEndian>(options.block_size)?;
    file.write_u64::<NativeEndian>(section_count)?;
    match &nodes {
        RTreeChildren::DataSections(sections) => {
            file.write_u32::<NativeEndian>(sections.first().unwrap().chrom)?;
            file.write_u32::<NativeEndian>(sections.first().unwrap().start)?;
            file.write_u32::<NativeEndian>(sections.last().unwrap().chrom)?;
            file.write_u32::<NativeEndian>(sections.last().unwrap().end)?;
        }
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

    let mut next_offset = file.tell()?;
    for level in (0..=levels).rev() {
        if level > 0 {
            next_offset += index_offsets[level - 1];
        }
        write_tree(file, &nodes, levels, level, next_offset, options)?;
    }

    Ok(())
}

pub(crate) fn write_zooms(
    mut file: &mut BufWriter<File>,
    zooms: Vec<ZoomInfo>,
    data_size: u64,
    options: BBIWriteOptions,
) -> io::Result<Vec<ZoomHeader>> {
    let mut zoom_entries: Vec<ZoomHeader> = vec![];
    let mut zoom_count = 0;
    let mut last_zoom_section_count = u64::max_value();
    for zoom in zooms {
        let zoom_file = zoom.data;
        let zoom_size = zoom_file.len()?;
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

        zoom_file.expect_closed_write(&mut file)?;
        let zoom_index_offset = file.tell()?;
        //println!("Zoom {:?}, data: {:?}, offset {:?}", zoom.resolution, zoom_data_offset, zoom_index_offset);
        assert_eq!(zoom_index_offset - zoom_data_offset, zoom_size);
        write_rtreeindex(&mut file, nodes, levels, total_sections, options)?;

        zoom_entries.push(ZoomHeader {
            reduction_level: zoom.resolution,
            data_offset: zoom_data_offset,
            index_offset: zoom_index_offset,
            index_tree_offset: None,
        });

        zoom_count += 1;
        if zoom_count >= options.max_zooms {
            break;
        }
    }

    Ok(zoom_entries)
}

/// Potential states encountered when reading `ChromData`
pub enum ChromDataState<ChromOutput, Error> {
    /// We've encountered a new chromosome
    NewChrom(ChromOutput),
    Finished,
    Error(Error),
}

/// An opaque key to indicate an processing chromosome
pub struct ChromProcessingKey(pub(crate) u32);

/// Effectively like an Iterator of chromosome data
pub trait ChromData: Sized {
    type Values: ChromValues;

    fn advance<
        State,
        F: FnMut(
            String,
            Self::Values,
            &mut State,
        ) -> Result<
            ChromProcessingKey,
            ProcessChromError<<Self::Values as ChromValues>::Error>,
        >,
    >(
        &mut self,
        do_read: &mut F,
        state: &mut State,
    ) -> Result<
        ChromDataState<ChromProcessingKey, <Self::Values as ChromValues>::Error>,
        ProcessChromError<<Self::Values as ChromValues>::Error>,
    >;
}

// Zooms have to be double-buffered: first because chroms could be processed in parallel and second because we don't know the offset of each zoom immediately
type ZoomValue = (
    Vec<crossbeam_channel::IntoIter<Section>>,
    TempFileBuffer<File>,
    Option<TempFileBufferWriter<File>>,
);
type Data<Error> = (
    crossbeam_channel::Receiver<Section>,
    TempFileBuffer<BufWriter<File>>,
    futures::future::RemoteHandle<Result<(usize, usize), ProcessChromError<Error>>>,
    Vec<TempZoomInfo<Error>>,
);
type DataWithoutzooms<Error> = (
    crossbeam_channel::Receiver<Section>,
    TempFileBuffer<BufWriter<File>>,
    futures::future::RemoteHandle<Result<(usize, usize), ProcessChromError<Error>>>,
);

async fn write_chroms_with_zooms<Err: Error + Send + 'static>(
    mut file: BufWriter<File>,
    mut zooms_map: BTreeMap<u32, ZoomValue>,
    mut receiver: futures_mpsc::UnboundedReceiver<Data<Err>>,
) -> Result<
    (
        BufWriter<File>,
        usize,
        Vec<crossbeam_channel::IntoIter<Section>>,
        BTreeMap<u32, ZoomValue>,
    ),
    ProcessChromError<Err>,
> {
    let mut section_iter = vec![];
    let mut max_uncompressed_buf_size = 0;
    loop {
        let read = receiver.next().await;
        let (sections, mut data, data_write_future, mut zooms) = match read {
            None => break,
            Some(read) => read,
        };
        // If we concurrently processing multiple chromosomes, the section buffer might have written some or all to a separate file
        // Switch that processing output to the real file
        data.switch(file);
        for TempZoomInfo {
            resolution: size,
            data: buf,
            ..
        } in zooms.iter_mut()
        {
            let zoom = zooms_map.get_mut(size).unwrap();
            let writer = zoom.2.take().unwrap();
            buf.switch(writer);
        }

        // All the futures are actually just handles, so these are purely for the result
        let (_num_sections, uncompressed_buf_size) = data_write_future.await?;
        max_uncompressed_buf_size = max_uncompressed_buf_size.max(uncompressed_buf_size);
        section_iter.push(sections.into_iter());
        file = data.await_real_file();

        for TempZoomInfo {
            resolution,
            data_write_future,
            data,
            sections,
        } in zooms.into_iter()
        {
            let zoom = zooms_map.get_mut(&resolution).unwrap();
            let data_write_data = data_write_future.await;
            let (_num_sections, uncompressed_buf_size) = match data_write_data {
                Ok(d) => d,
                Err(e) => return Err(e),
            };
            max_uncompressed_buf_size = max_uncompressed_buf_size.max(uncompressed_buf_size);
            zoom.0.push(sections.into_iter());
            zoom.2.replace(data.await_real_file());
        }
    }

    Ok((file, max_uncompressed_buf_size, section_iter, zooms_map))
}

async fn write_chroms_without_zooms<Err: Error + Send + 'static>(
    mut file: BufWriter<File>,
    mut receiver: futures_mpsc::UnboundedReceiver<DataWithoutzooms<Err>>,
) -> Result<
    (
        BufWriter<File>,
        usize,
        Vec<crossbeam_channel::IntoIter<Section>>,
    ),
    ProcessChromError<Err>,
> {
    let mut section_iter = vec![];
    let mut max_uncompressed_buf_size = 0;
    loop {
        let read = receiver.next().await;
        let (sections, mut data, data_write_future) = match read {
            None => break,
            Some(read) => read,
        };
        // If we concurrently processing multiple chromosomes, the section buffer might have written some or all to a separate file
        // Switch that processing output to the real file
        data.switch(file);

        // All the futures are actually just handles, so these are purely for the result
        let (_num_sections, uncompressed_buf_size) = data_write_future.await?;
        max_uncompressed_buf_size = max_uncompressed_buf_size.max(uncompressed_buf_size);
        section_iter.push(sections.into_iter());
        file = data.await_real_file();
    }

    Ok((file, max_uncompressed_buf_size, section_iter))
}

pub(crate) fn write_vals<
    Values: ChromValues,
    V: ChromData<Values = Values>,
    Fut: Future<Output = Result<Summary, ProcessChromError<Values::Error>>>,
    G: Fn(
        Vec<(u32, ChromProcessingInputSectionChannel)>,
        ChromProcessingInputSectionChannel,
        u32,
        BBIWriteOptions,
        Handle,
        Values,
        String,
        u32,
    ) -> Fut,
>(
    mut vals_iter: V,
    file: BufWriter<File>,
    options: BBIWriteOptions,
    process_chrom: G,
    runtime: Runtime,
    chrom_sizes: HashMap<String, u32>,
) -> Result<
    (
        IdMap,
        Summary,
        BufWriter<File>,
        Flatten<vec::IntoIter<crossbeam_channel::IntoIter<Section>>>,
        Vec<ZoomInfo>,
        usize,
    ),
    ProcessChromError<Values::Error>,
> {
    let zooms_map: BTreeMap<u32, ZoomValue> =
        std::iter::successors(Some(options.initial_zoom_size), |z| Some(z * 4))
            .take(options.max_zooms as usize)
            .map(|size| {
                let section_iter = vec![];
                let (buf, write): (TempFileBuffer<File>, TempFileBufferWriter<File>) =
                    TempFileBuffer::new(options.inmemory);
                let value = (section_iter, buf, Some(write));
                (size, value)
            })
            .collect();

    let mut chrom_ids = IdMap::default();

    let mut key = 0;
    let mut output: BTreeMap<u32, _> = BTreeMap::new();

    let mut summary: Option<Summary> = None;
    let (send, recv) = futures_mpsc::unbounded();
    let write_fut = write_chroms_with_zooms(file, zooms_map, recv);

    let setup_chrom = || {
        let (ftx, sections_handle, buf, section_receiver) =
            future_channel(options.channel_size, runtime.handle(), options.inmemory);

        let (zoom_infos, zooms_channels) = {
            let mut zoom_infos = Vec::with_capacity(options.max_zooms as usize);
            let mut zooms_channels = Vec::with_capacity(options.max_zooms as usize);

            let zoom_sizes =
                std::iter::successors(Some(options.initial_zoom_size), |z| Some(z * 4))
                    .take(options.max_zooms as usize);
            for size in zoom_sizes {
                let (ftx, handle, buf, section_receiver) =
                    future_channel(options.channel_size, runtime.handle(), options.inmemory);
                let zoom_info = TempZoomInfo {
                    resolution: size,
                    data_write_future: Box::new(handle),
                    data: buf,
                    sections: section_receiver,
                };
                zoom_infos.push(zoom_info);
                zooms_channels.push((size, ftx));
            }
            (zoom_infos, zooms_channels)
        };

        match send.unbounded_send((section_receiver, buf, sections_handle, zoom_infos)) {
            Ok(_) => {}
            Err(_) => panic!("Expected to always send."),
        }

        (zooms_channels, ftx)
    };
    let mut do_read = |chrom: String,
                       data: _,
                       output: &mut BTreeMap<u32, _>|
     -> Result<ChromProcessingKey, ProcessChromError<_>> {
        let length = match chrom_sizes.get(&chrom) {
            Some(length) => *length,
            None => {
                return Err(ProcessChromError::InvalidChromosome(format!(
                    "Input bedGraph contains chromosome that isn't in the input chrom sizes: {}",
                    chrom
                )));
            }
        };
        // Make a new id for the chromosome
        let chrom_id = chrom_ids.get_id(&chrom);

        let (zooms_channels, ftx) = setup_chrom();

        let fut = process_chrom(
            zooms_channels,
            ftx,
            chrom_id,
            options,
            runtime.handle().clone(),
            data,
            chrom,
            length,
        );

        let curr_key = key;
        key += 1;

        output.insert(curr_key, fut);

        Ok(ChromProcessingKey(curr_key))
    };

    let (write_fut, write_fut_handle) = write_fut.remote_handle();
    runtime.spawn(write_fut);
    loop {
        match vals_iter.advance(&mut do_read, &mut output)? {
            ChromDataState::NewChrom(read) => {
                let fut = output.remove(&read.0).unwrap();
                let chrom_summary = runtime.block_on(fut)?;
                match &mut summary {
                    None => summary = Some(chrom_summary),
                    Some(summary) => {
                        summary.total_items += chrom_summary.total_items;
                        summary.bases_covered += chrom_summary.bases_covered;
                        summary.min_val = summary.min_val.min(chrom_summary.min_val);
                        summary.max_val = summary.max_val.max(chrom_summary.max_val);
                        summary.sum += chrom_summary.sum;
                        summary.sum_squares += chrom_summary.sum_squares;
                    }
                }
            }
            ChromDataState::Finished => break,
            ChromDataState::Error(err) => return Err(ProcessChromError::SourceError(err)),
        }
    }
    drop(send);

    let summary_complete = summary.unwrap_or(Summary {
        total_items: 0,
        bases_covered: 0,
        min_val: 0.0,
        max_val: 0.0,
        sum: 0.0,
        sum_squares: 0.0,
    });

    let (file, max_uncompressed_buf_size, section_iter, zooms_map) =
        runtime.block_on(write_fut_handle)?;

    let zoom_infos: Vec<ZoomInfo> = zooms_map
        .into_iter()
        .map(|(size, zoom)| {
            drop(zoom.2);
            let sections = zoom.0.into_iter().flatten();
            ZoomInfo {
                resolution: size,
                data: zoom.1,
                sections,
            }
        })
        .collect();
    let section_iter = section_iter.into_iter().flatten();
    Ok((
        chrom_ids,
        summary_complete,
        file,
        section_iter,
        zoom_infos,
        max_uncompressed_buf_size,
    ))
}

pub(crate) fn write_vals_no_zoom<
    Values: ChromValues,
    V: ChromData<Values = Values>,
    Fut: Future<Output = Result<(Summary, Vec<(u64, u64)>), ProcessChromError<Values::Error>>>
        + Send
        + 'static,
    G: Fn(
        ChromProcessingInputSectionChannel,
        u32,
        BBIWriteOptions,
        Handle,
        Values,
        String,
        u32,
    ) -> Fut,
>(
    mut vals_iter: V,
    file: BufWriter<File>,
    options: BBIWriteOptions,
    process_chrom: G,
    runtime: &Runtime,
    chrom_sizes: HashMap<String, u32>,
) -> Result<
    (
        IdMap,
        Summary,
        BTreeMap<u64, u64>,
        BufWriter<File>,
        Flatten<vec::IntoIter<crossbeam_channel::IntoIter<Section>>>,
        usize,
    ),
    ProcessChromError<Values::Error>,
> {
    let total_zoom_counts = std::iter::successors(Some(10), |z: &u64| Some((*z).saturating_mul(4)))
        .take_while(|z| *z < u64::MAX)
        .map(|z| (z, 0));
    let mut total_zoom_counts: BTreeMap<u64, u64> = BTreeMap::from_iter(total_zoom_counts);

    let mut chrom_ids = IdMap::default();

    let mut key = 0;
    let mut output: BTreeMap<u32, _> = BTreeMap::new();

    let mut summary: Option<Summary> = None;
    let (send, recv) = futures_mpsc::unbounded();
    let write_fut = write_chroms_without_zooms::<Values::Error>(file, recv);

    let setup_chrom = || {
        let (ftx, sections_handle, buf, section_receiver) =
            future_channel(options.channel_size, runtime.handle(), options.inmemory);

        match send.unbounded_send((section_receiver, buf, sections_handle)) {
            Ok(_) => {}
            Err(_) => panic!("Expected to always send."),
        }

        ftx
    };
    let mut do_read = |chrom: String,
                       data: _,
                       output: &mut BTreeMap<u32, _>|
     -> Result<ChromProcessingKey, ProcessChromError<_>> {
        let length = match chrom_sizes.get(&chrom) {
            Some(length) => *length,
            None => {
                return Err(ProcessChromError::InvalidChromosome(format!(
                    "Input bedGraph contains chromosome that isn't in the input chrom sizes: {}",
                    chrom
                )));
            }
        };
        // Make a new id for the chromosome
        let chrom_id = chrom_ids.get_id(&chrom);

        let ftx = setup_chrom();

        let fut = process_chrom(
            ftx,
            chrom_id,
            options,
            runtime.handle().clone(),
            data,
            chrom,
            length,
        );

        let curr_key = key;
        key += 1;

        output.insert(curr_key, fut);

        Ok(ChromProcessingKey(curr_key))
    };

    let (write_fut, write_fut_handle) = write_fut.remote_handle();
    runtime.spawn(write_fut);
    loop {
        match vals_iter.advance(&mut do_read, &mut output)? {
            ChromDataState::NewChrom(read) => {
                let fut = output.remove(&read.0).unwrap();
                let (chrom_summary, zoom_counts) = runtime.block_on(fut)?;

                match &mut summary {
                    None => summary = Some(chrom_summary),
                    Some(summary) => {
                        summary.total_items += chrom_summary.total_items;
                        summary.bases_covered += chrom_summary.bases_covered;
                        summary.min_val = summary.min_val.min(chrom_summary.min_val);
                        summary.max_val = summary.max_val.max(chrom_summary.max_val);
                        summary.sum += chrom_summary.sum;
                        summary.sum_squares += chrom_summary.sum_squares;
                    }
                }

                let zoom_count_map = BTreeMap::from_iter(zoom_counts.into_iter());
                for zoom_count in total_zoom_counts.iter_mut() {
                    let chrom_zoom_count = zoom_count_map.get(&zoom_count.0).copied().unwrap_or(1);
                    *zoom_count.1 += chrom_zoom_count;
                }
            }
            ChromDataState::Finished => break,
            ChromDataState::Error(err) => return Err(ProcessChromError::SourceError(err)),
        }
    }
    drop(send);

    let summary_complete = summary.unwrap_or(Summary {
        total_items: 0,
        bases_covered: 0,
        min_val: 0.0,
        max_val: 0.0,
        sum: 0.0,
        sum_squares: 0.0,
    });

    let (file, max_uncompressed_buf_size, section_iter) = runtime.block_on(write_fut_handle)?;

    let section_iter = section_iter.into_iter().flatten();
    Ok((
        chrom_ids,
        summary_complete,
        total_zoom_counts,
        file,
        section_iter,
        max_uncompressed_buf_size,
    ))
}

pub(crate) fn write_zoom_vals<
    Values: ChromValues,
    V: ChromData<Values = Values>,
    Fut: Future<Output = Result<(), ProcessChromError<Values::Error>>> + Send + 'static,
    G: Fn(Vec<(u32, ChromProcessingInputSectionChannel)>, u32, BBIWriteOptions, Handle, Values) -> Fut,
>(
    mut vals_iter: V,
    options: BBIWriteOptions,
    process_chrom_zoom: G,
    runtime: &Runtime,
    chrom_ids: &HashMap<String, u32>,
    average_size: u32,
    zoom_counts: BTreeMap<u64, u64>,
    mut file: BufWriter<File>,
    data_size: u64,
) -> Result<(BufWriter<File>, Vec<ZoomHeader>, usize), ProcessChromError<Values::Error>> {
    // Zooms have to be double-buffered: first because chroms could be processed in parallel and second because we don't know the offset of each zoom immediately
    type ZoomValue = (
        Vec<crossbeam_channel::IntoIter<Section>>,
        TempFileBuffer<BufWriter<File>>,
        Option<TempFileBufferWriter<BufWriter<File>>>,
    );

    pub(crate) struct TempZoomInfo<SourceError: Error> {
        pub resolution: u32,
        pub data_write_future: Box<
            dyn Future<Output = Result<(usize, usize), ProcessChromError<SourceError>>>
                + Send
                + Unpin,
        >,
        pub data: TempFileBuffer<TempFileBufferWriter<BufWriter<File>>>,
        pub sections: crossbeam_channel::Receiver<Section>,
    }

    let min_first_zoom_size = average_size.max(10) * 4;
    let mut zooms_map: BTreeMap<u32, ZoomValue> = zoom_counts
        .into_iter()
        .skip_while(|z| z.0 > min_first_zoom_size as u64)
        .skip_while(|z| {
            let mut reduced_size = z.1 * 32;
            if options.compress {
                reduced_size /= 2; // Estimate as kent does
            }
            reduced_size as u64 > data_size / 2
        })
        .take(options.max_zooms as usize)
        .map(|size| {
            let section_iter = vec![];
            let (buf, write) = TempFileBuffer::new(options.inmemory);
            let value = (section_iter, buf, Some(write));
            (size.0 as u32, value)
        })
        .collect();
    let resolutions: Vec<_> = zooms_map.keys().copied().collect();

    let first_zoom_data_offset = file.tell()?;
    // We can immediately start to write to the file the first zoom
    match zooms_map.first_entry() {
        Some(mut first) => first.get_mut().1.switch(file),
        None => return Ok((file, vec![], 0)),
    }

    let mut max_uncompressed_buf_size = 0;

    let mut key = 0;
    let mut output: BTreeMap<u32, _> = BTreeMap::new();

    let mut do_read = |chrom: String,
                       data: _,
                       output: &mut BTreeMap<u32, _>|
     -> Result<ChromProcessingKey, ProcessChromError<_>> {
        // Make a new id for the chromosome
        let chrom_id = *chrom_ids
            .get(&chrom)
            .expect("Should not have seen a new chrom.");

        let (zoom_infos, zooms_channels) = {
            let mut zoom_infos = Vec::with_capacity(options.max_zooms as usize);
            let mut zooms_channels = Vec::with_capacity(options.max_zooms as usize);

            for size in resolutions.iter().copied() {
                let (ftx, handle, buf, section_receiver) =
                    future_channel(options.channel_size, runtime.handle(), options.inmemory);
                let zoom_info = TempZoomInfo {
                    resolution: size,
                    data_write_future: Box::new(handle),
                    data: buf,
                    sections: section_receiver,
                };
                zoom_infos.push(zoom_info);
                zooms_channels.push((size, ftx));
            }
            (zoom_infos, zooms_channels)
        };

        let (f_remote, f_handle) = process_chrom_zoom(
            zooms_channels,
            chrom_id,
            options,
            runtime.handle().clone(),
            data,
        )
        .remote_handle();
        runtime.spawn(f_remote);

        let curr_key = key;
        key += 1;

        output.insert(curr_key, (f_handle, zoom_infos));

        Ok(ChromProcessingKey(curr_key))
    };

    loop {
        match vals_iter.advance(&mut do_read, &mut output)? {
            ChromDataState::NewChrom(read) => {
                let read = output.remove(&read.0).unwrap();
                let (process_future, mut zooms) = read;
                // For each zoom, switch the current chromosome to write to the actual zoom file
                for TempZoomInfo {
                    resolution: size,
                    data,
                    ..
                } in zooms.iter_mut()
                {
                    let zoom = zooms_map.get_mut(size).unwrap();
                    let writer = zoom.2.take().unwrap();
                    data.switch(writer);
                }

                runtime.block_on(process_future)?;

                for TempZoomInfo {
                    resolution,
                    data_write_future,
                    data,
                    sections,
                } in zooms.into_iter()
                {
                    // First, we need to make sure that all the sections that were queued to encode have been written
                    let data_write_data = runtime.block_on(data_write_future);
                    let (_num_sections, uncompressed_buf_size) = match data_write_data {
                        Ok(d) => d,
                        Err(e) => return Err(e),
                    };
                    max_uncompressed_buf_size =
                        max_uncompressed_buf_size.max(uncompressed_buf_size);

                    let zoom = zooms_map.get_mut(&resolution).unwrap();
                    // Add the section data to the zoom
                    zoom.0.push(sections.into_iter());
                    // Replace the zoom file again
                    zoom.2.replace(data.await_real_file());
                }
            }
            ChromDataState::Finished => break,
            ChromDataState::Error(err) => return Err(ProcessChromError::SourceError(err)),
        }
    }

    let mut zoom_entries = Vec::with_capacity(zooms_map.len());
    let mut zooms_map_iter = zooms_map.into_iter();

    // Since the first zoom has already been written to the file, no need to
    let first_zoom = zooms_map_iter
        .next()
        .expect("Should have at least one zoom");
    // First, we can drop the writer - no more data
    drop(first_zoom.1 .2);
    let first_zoom_sections = first_zoom.1 .0.into_iter().flatten();
    let mut current_offset = first_zoom_data_offset;
    let sections_iter = first_zoom_sections.map(|mut section| {
        // TODO: assumes contiguous, see note for primary data
        section.offset = current_offset;
        current_offset += section.size;
        section
    });
    // First zoom has already switched, real data
    file = first_zoom.1 .1.await_real_file();
    // Generate the rtree index
    let (nodes, levels, total_sections) = get_rtreeindex(sections_iter, options);
    let first_zoom_index_offset = file.tell()?;
    write_rtreeindex(&mut file, nodes, levels, total_sections, options)?;
    zoom_entries.push(ZoomHeader {
        reduction_level: first_zoom.0,
        data_offset: first_zoom_data_offset,
        index_offset: first_zoom_index_offset,
        index_tree_offset: None,
    });

    while let Some(mut zoom) = zooms_map_iter.next() {
        let zoom_data_offset = file.tell()?;
        // First, we can drop the writer - no more data
        drop(zoom.1 .2);
        let zoom_sections = zoom.1 .0.into_iter().flatten();
        let mut current_offset = zoom_data_offset;
        let sections_iter = zoom_sections.map(|mut section| {
            // TODO: assumes contiguous, see note for primary data
            section.offset = current_offset;
            current_offset += section.size;
            section
        });
        // Subsequence zooms have not switched to real file
        zoom.1 .1.switch(file);
        file = zoom.1 .1.await_real_file();
        // Generate the rtree index
        let (nodes, levels, total_sections) = get_rtreeindex(sections_iter, options);
        let zoom_index_offset = file.tell()?;
        write_rtreeindex(&mut file, nodes, levels, total_sections, options)?;
        zoom_entries.push(ZoomHeader {
            reduction_level: first_zoom.0,
            data_offset: zoom_data_offset,
            index_offset: zoom_index_offset,
            index_tree_offset: None,
        });
    }

    Ok((file, zoom_entries, max_uncompressed_buf_size))
}

pub(crate) fn write_mid<E: Error>(
    file: &mut BufWriter<File>,
    pre_data: u64,
    raw_sections_iter: impl Iterator<Item = Section>,
    chrom_sizes: HashMap<String, u32>,
    chrom_ids: &HashMap<String, u32>,
    options: BBIWriteOptions,
) -> Result<(u64, u64, u64, u64), ProcessChromError<E>> {
    let data_size = file.tell()? - pre_data;
    let mut current_offset = pre_data;
    let sections_iter = raw_sections_iter.map(|mut section| {
        // TODO: this assumes that all the data is contiguous
        // This will fail if we ever space the sections in any way
        section.offset = current_offset;
        current_offset += section.size;
        section
    });

    // This deviates slighly from the layout of bigBeds generated from kent tools (but are 100%)
    // compatible. In kent tools, the chrom tree is written *before* the data.
    // However, in order to do this, the data must be read *twice*, which we don't want. Luckily,
    // the chrom tree offset is stored to be seeked to before read, so
    // it doesn't matter where in the file these are placed. The one caveat to this is in any
    // caching (most commonly over-the-network): if the caching "chunk" size is large enough to
    // cover either/both the autosql and chrom tree with the start of the file, then this may
    // cause two separate "chunks" (and therefore, e.g., network requests), compared to one.

    // Since the chrom tree is read before the index, we put this before the full data index
    // Therefore, there is a higher likelihood that the udc file will only need one read for
    // chrom tree + full data index.
    let chrom_index_start = file.tell()?;
    write_chrom_tree(file, chrom_sizes, &chrom_ids)?;

    let index_start = file.tell()?;
    let (nodes, levels, total_sections) = get_rtreeindex(sections_iter, options);
    write_rtreeindex(file, nodes, levels, total_sections, options)?;

    Ok((data_size, chrom_index_start, index_start, total_sections))
}

async fn write_data<W: Write, SourceError: Error + Send>(
    mut data_file: W,
    section_sender: crossbeam_channel::Sender<Section>,
    mut frx: futures_mpsc::Receiver<impl Future<Output = io::Result<(SectionData, usize)>> + Send>,
) -> Result<(usize, usize), ProcessChromError<SourceError>> {
    let mut current_offset = 0;
    let mut total = 0;
    let mut max_uncompressed_buf_size = 0;
    while let Some(section_raw) = frx.next().await {
        let (section, uncompressed_buf_size): (SectionData, usize) = section_raw.await?;
        max_uncompressed_buf_size = max_uncompressed_buf_size.max(uncompressed_buf_size);
        total += 1;
        let size = section.data.len() as u64;
        data_file.write_all(&section.data)?;
        section_sender
            .send(Section {
                chrom: section.chrom,
                start: section.start,
                end: section.end,
                offset: current_offset,
                size,
            })
            .expect("Couldn't send section.");
        current_offset += size;
    }
    Ok((total, max_uncompressed_buf_size))
}

pub(crate) fn future_channel<Err: Error + Send + 'static, R: Write + Send + 'static>(
    channel_size: usize,
    runtime: &Handle,
    inmemory: bool,
) -> (
    futures_mpsc::Sender<
        Pin<Box<dyn Future<Output = Result<(SectionData, usize), io::Error>> + Send>>,
    >,
    futures::future::RemoteHandle<Result<(usize, usize), ProcessChromError<Err>>>,
    TempFileBuffer<R>,
    crossbeam_channel::Receiver<Section>,
) {
    let (ftx, frx) = channel(channel_size);
    let (buf, write) = TempFileBuffer::new(inmemory);
    let file = BufWriter::new(write);

    let (section_sender, section_receiver) = unbounded();
    let (sections_remote, sections_handle) = write_data(file, section_sender, frx).remote_handle();
    runtime.spawn(sections_remote);
    (ftx, sections_handle, buf, section_receiver)
}

#[cfg(all(test, feature = "read"))]
mod tests {
    use byteordered::Endianness;

    use crate::{read_cir_tree_header, search_cir_tree_inner};

    use super::*;
    use std::io::{BufReader, Cursor};

    #[test]
    fn test_rtreeindex() -> io::Result<()> {
        const MAX_BASES: u32 = 256 * 256 * 256;
        let mut chrom = 0;
        let mut start = 0;
        let iter = std::iter::from_fn(|| {
            let curr_chrom = chrom;
            let curr_start = start;
            start += 1;
            if start >= MAX_BASES {
                start = 0;
                chrom += 1;
            }
            Some(Section {
                chrom: curr_chrom,
                start: curr_start,
                end: curr_start + 1,
                offset: curr_start as u64,
                size: 1,
            })
        });
        let mut options = BBIWriteOptions::default();
        options.block_size = 5;
        let (tree, levels, total_sections) = get_rtreeindex(iter.take(126), options);

        let mut data = Vec::<u8>::new();
        let mut cursor = Cursor::new(&mut data);
        let mut bufwriter = BufWriter::new(&mut cursor);
        write_rtreeindex(&mut bufwriter, tree, levels, total_sections, options)?;

        drop(bufwriter);
        drop(cursor);

        let mut cursor = Cursor::new(&mut data);
        let mut file = BufReader::new(&mut cursor);

        read_cir_tree_header(Endianness::native(), &mut file).unwrap();

        let blocks =
            search_cir_tree_inner(Endianness::native(), &mut file, 48, 0, 0, MAX_BASES).unwrap();

        let mut chrom = 0;
        let mut start = 0;
        let iter = std::iter::from_fn(|| {
            let curr_chrom = chrom;
            let curr_start = start;
            start += 1;
            if start >= MAX_BASES {
                start = 0;
                chrom += 1;
            }
            Some(Section {
                chrom: curr_chrom,
                start: curr_start,
                end: curr_start + 1,
                offset: curr_start as u64,
                size: 1,
            })
        });
        iter.zip(blocks.into_iter())
            .for_each(|(a, b)| assert_eq!(a.offset, b.offset));
        Ok(())
    }
}
