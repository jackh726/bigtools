use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::pin::Pin;

use byteorder::{NativeEndian, WriteBytesExt};
use thiserror::Error;

use futures::channel::mpsc::{channel, Receiver};
use futures::executor::ThreadPool;
use futures::future::{Future, FutureExt};
use futures::stream::StreamExt;
use futures::task::SpawnExt;
use futures::try_join;

use serde::{Deserialize, Serialize};

use crate::utils::chromvalues::ChromValues;
use crate::utils::filebufferedchannel;
use crate::utils::idmap::IdMap;
use crate::utils::tell::Tell;
use crate::utils::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};

use crate::bbi::{Summary, ZoomHeader, ZoomRecord, CHROM_TREE_MAGIC, CIR_TREE_MAGIC};

pub(crate) struct ZoomInfo {
    resolution: u32,
    data: File,
    sections: Box<dyn Iterator<Item = Section>>,
}

#[derive(Debug)]
pub(crate) struct SectionData {
    pub(crate) chrom: u32,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) data: Vec<u8>,
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Section {
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

#[derive(Copy, Clone)]
pub enum InputSortType {
    ALL,
    START,
    // TODO
    //NONE,
}

#[derive(Copy, Clone)]
pub struct BBIWriteOptions {
    pub compress: bool,
    pub items_per_slot: u32,
    pub block_size: u32,
    pub initial_zoom_size: u32,
    pub max_zooms: u32,
    pub input_sort_type: InputSortType,
    pub channel_size: usize,
}

impl Default for BBIWriteOptions {
    fn default() -> Self {
        BBIWriteOptions {
            compress: true,
            items_per_slot: 1024,
            block_size: 256,
            initial_zoom_size: 160,
            max_zooms: 10,
            input_sort_type: InputSortType::ALL,
            channel_size: 100,
        }
    }
}

#[derive(Error, Debug)]
pub enum ProcessChromError<SourceError> {
    #[error("{}", .0)]
    InvalidInput(String),
    #[error("{}", .0)]
    InvalidChromosome(String),
    #[error("{}", .0)]
    IoError(#[from] io::Error),
    #[error("SourceError")]
    SourceError(SourceError),
}

impl<W, S> From<io::IntoInnerError<W>> for ProcessChromError<S> {
    fn from(error: io::IntoInnerError<W>) -> Self {
        ProcessChromError::IoError(error.into())
    }
}

pub struct TempZoomInfo<SourceError> {
    pub resolution: u32,
    pub data_write_future: Box<
        dyn Future<Output = Result<(usize, usize), ProcessChromError<SourceError>>> + Send + Unpin,
    >,
    pub data: TempFileBuffer<TempFileBufferWriter<File>>,
    pub sections: filebufferedchannel::Receiver<Section>,
}

pub(crate) type ChromProcessingInputSectionChannel = futures::channel::mpsc::Sender<
    Pin<Box<dyn Future<Output = io::Result<(SectionData, usize)>> + Send>>,
>;
pub(crate) struct ChromProcessingInput {
    pub(crate) zooms_channels: Vec<ChromProcessingInputSectionChannel>,
    pub(crate) ftx: ChromProcessingInputSectionChannel,
}

pub struct ChromProcessingOutput<SourceError> {
    pub sections: filebufferedchannel::Receiver<Section>,
    pub data: TempFileBuffer<File>,
    pub data_write_future: Box<
        dyn Future<Output = Result<(usize, usize), ProcessChromError<SourceError>>> + Send + Unpin,
    >,
    pub zooms: Vec<TempZoomInfo<SourceError>>,
}

pub type WriteSummaryFuture<SourceError> =
    Pin<Box<dyn Future<Output = Result<Summary, ProcessChromError<SourceError>>> + Send>>;

const MAX_ZOOM_LEVELS: usize = 10;

pub(crate) fn write_blank_headers(file: &mut BufWriter<File>) -> io::Result<()> {
    file.seek(SeekFrom::Start(0))?;
    // Common header
    file.write_all(&[0; 64])?;
    // Zoom levels
    file.write_all(&[0; MAX_ZOOM_LEVELS * 24])?;

    Ok(())
}

pub(crate) fn write_chrom_tree(
    file: &mut BufWriter<File>,
    chrom_sizes: std::collections::HashMap<String, u32>,
    chrom_ids: &std::collections::HashMap<String, u32>,
) -> io::Result<()> {
    let mut chroms: Vec<&String> = chrom_ids.keys().collect();
    chroms.sort();
    //println!("Used chroms {:?}", chroms);

    let item_count = chroms.len() as u64;
    // TODO: for now, always just use the length of chroms (if less than 256). This means we don't have to implement writing non-leaf nodes for now...
    // TODO: make this configurable
    let block_size = std::cmp::max(256, item_count) as u32;
    let max_bytes = chroms.iter().map(|a| a.len() as u32).fold(0, u32::max);

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
    for chrom in chroms {
        let key_bytes = &mut vec![0u8; max_bytes as usize];
        let chrom_bytes = chrom.as_bytes();
        key_bytes[..chrom_bytes.len()].copy_from_slice(chrom_bytes);
        file.write_all(key_bytes)?;
        let id = *chrom_ids
            .get(chrom)
            .expect("Internal error. (Chrom not found).");
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
        io::copy(&mut buf_reader, &mut file)?;
        let zoom_index_offset = file.tell()?;
        //println!("Zoom {:?}, data: {:?}, offset {:?}", zoom.resolution, zoom_data_offset, zoom_index_offset);
        assert_eq!(zoom_index_offset - zoom_data_offset, zoom_size);
        write_rtreeindex(&mut file, nodes, levels, total_sections, options)?;

        zoom_entries.push(ZoomHeader {
            reduction_level: zoom.resolution,
            data_offset: zoom_data_offset,
            index_offset: zoom_index_offset,
        });

        zoom_count += 1;
        if zoom_count >= options.max_zooms {
            break;
        }
    }

    Ok(zoom_entries)
}

/// Potential states encountered when reading `ChromData`
pub enum ChromDataState<Error> {
    /// We've encountered a new chromosome
    NewChrom(ChromProcessingFnOutput<Error>),
    Finished,
    Error(Error),
}

/// Effectively like an Iterator of chromosome data
pub trait ChromData<E: From<io::Error>>: Sized {
    type Output: ChromValues;
    fn advance<
        F: FnMut(
            String,
            Self::Output,
        ) -> Result<ChromProcessingFnOutput<<Self::Output as ChromValues>::Error>, E>,
    >(
        &mut self,
        do_read: &mut F,
    ) -> Result<ChromDataState<<Self::Output as ChromValues>::Error>, E>;
}

pub struct ChromProcessingFnOutput<Error>(pub(crate) WriteSummaryFuture<Error>, pub(crate) ChromProcessingOutput<Error>);

pub(crate) async fn write_vals<
    Values: ChromValues,
    V: ChromData<ProcessChromError<Values::Error>, Output = Values>,
    Fut: Future<Output = Result<Summary, ProcessChromError<Values::Error>>> + Send + 'static,
    G: Fn(ChromProcessingInput, u32, BBIWriteOptions, ThreadPool, Values, String, u32) -> Fut,
>(
    mut vals_iter: V,
    file: BufWriter<File>,
    options: BBIWriteOptions,
    process_chrom: G,
    mut pool: ThreadPool,
    chrom_sizes: HashMap<String, u32>,
) -> Result<
    (
        IdMap,
        Summary,
        BufWriter<File>,
        Box<dyn Iterator<Item = Section> + 'static>,
        Vec<ZoomInfo>,
        usize,
    ),
    ProcessChromError<Values::Error>,
> {
    // Zooms have to be double-buffered: first because chroms could be processed in parallel and second because we don't know the offset of each zoom immediately
    type ZoomValue = (
        Vec<Box<dyn Iterator<Item = Section>>>,
        TempFileBuffer<File>,
        Option<TempFileBufferWriter<File>>,
    );

    let mut zooms_map: BTreeMap<u32, ZoomValue> =
        std::iter::successors(Some(options.initial_zoom_size), |z| Some(z * 4))
            .take(options.max_zooms as usize)
            .map(|size| -> io::Result<_> {
                let section_iter: Vec<Box<dyn Iterator<Item = Section>>> = vec![];
                let (buf, write): (TempFileBuffer<File>, TempFileBufferWriter<File>) =
                    TempFileBuffer::new()?;
                let value = (section_iter, buf, Some(write));
                Ok((size, value))
            })
            .collect::<io::Result<_>>()?;

    let mut section_iter: Vec<filebufferedchannel::IntoIter<Section>> = vec![];
    let mut raw_file = file.into_inner()?;

    let mut summary: Option<Summary> = None;
    let mut max_uncompressed_buf_size = 0;

    let mut chrom_ids = IdMap::default();

    let mut do_read = |chrom: String,
                       data: _|
     -> Result<
        ChromProcessingFnOutput<<Values as ChromValues>::Error>,
        ProcessChromError<_>,
    > {
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

        // This converts a ChromValues (streaming iterator) to a (WriteSummaryFuture, ChromProcessingOutput).
        // This is a separate function so this can techincally be run for mulitple chromosomes simulatenously.
        // This is heavily multi-threaded using Futures. A brief summary:
        // - All reading from the ChromValues is done in a single future (process_chrom). This futures in charge of keeping track of sections (and zoom data).
        //   When a section is full, a Future is created to byte-encode the data and compress it (if compression is on). The same is true for zoom sections.
        //   This is the most CPU-intensive part of the entire write.
        // - The section futures are sent (in order) by channel to a separate future for the sole purpose of writing the (maybe compressed) section data to a `TempFileBuffer`.
        //   The data is written to a temporary file, since this may be happening in parallel (where we don't know the real file offset of these sections).
        //   `TempFileBuffer` allows "switching" to the real file once it's available (on the read side), so if the real file is available, I/O ops are not duplicated.
        //   Once the section data is written to the file, the file offset data is stored in a `FileBufferedChannel`. This is needed for the index.
        //   All of this is done for zoom sections too.
        //
        // The futures that are returned are only handles to remote futures that are spawned immediately on `pool`.
        let (procesing_input, processing_output) = setup_channels(&mut pool, options)?;

        let (f_remote, f_handle) = process_chrom(
            procesing_input,
            chrom_id,
            options,
            pool.clone(),
            data,
            chrom,
            length,
        )
        .remote_handle();
        pool.spawn(f_remote).expect("Couldn't spawn future.");
        Ok(ChromProcessingFnOutput(f_handle.boxed(), processing_output))
    };

    let chrom_ids = loop {
        match vals_iter.advance(&mut do_read)? {
            ChromDataState::NewChrom(read) => {
                let ChromProcessingFnOutput(
                    summary_future,
                    ChromProcessingOutput {
                        sections,
                        mut data,
                        data_write_future,
                        mut zooms,
                    },
                ) = read;
                // If we concurrently processing multiple chromosomes, the section buffer might have written some or all to a separate file
                // Switch that processing output to the real file
                data.switch(raw_file);
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
                let joined_future = match try_join!(summary_future, data_write_future) {
                    Ok(f) => f,
                    Err(e) => return Err(e),
                };
                let (chrom_summary, (_num_sections, uncompressed_buf_size)) = joined_future;
                max_uncompressed_buf_size = max_uncompressed_buf_size.max(uncompressed_buf_size);
                section_iter.push(sections.into_iter());
                raw_file = data.await_real_file();

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
                    max_uncompressed_buf_size =
                        max_uncompressed_buf_size.max(uncompressed_buf_size);
                    zoom.0.push(Box::new(sections.into_iter()));
                    zoom.2.replace(data.await_real_file());
                }

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
            ChromDataState::Finished => break chrom_ids,
            ChromDataState::Error(err) => return Err(ProcessChromError::SourceError(err)),
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

    let zoom_infos: Vec<ZoomInfo> = zooms_map
        .into_iter()
        .map(|(size, zoom)| {
            drop(zoom.2);
            let zoom_iter: Box<dyn Iterator<Item = Section> + 'static> =
                Box::new(zoom.0.into_iter().flatten());
            let closed_file = zoom.1.await_temp_file();
            ZoomInfo {
                resolution: size,
                data: closed_file,
                sections: zoom_iter,
            }
        })
        .collect();
    let section_iter: Box<dyn Iterator<Item = Section> + 'static> =
        Box::new(section_iter.into_iter().flatten());
    Ok((
        chrom_ids,
        summary_complete,
        BufWriter::new(raw_file),
        section_iter,
        zoom_infos,
        max_uncompressed_buf_size,
    ))
}

async fn write_data<W: Write, SourceError: Send>(
    mut data_file: W,
    mut section_sender: filebufferedchannel::Sender<Section>,
    mut frx: Receiver<impl Future<Output = io::Result<(SectionData, usize)>> + Send>,
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

/// Sets up the channels and write "threads" for the data and zoom sections
pub(crate) fn setup_channels<SourceError: Send + 'static>(
    pool: &mut ThreadPool,
    options: BBIWriteOptions,
) -> io::Result<(ChromProcessingInput, ChromProcessingOutput<SourceError>)> {
    let (ftx, frx) = channel(options.channel_size);

    let (sections_handle, buf, section_receiver) = {
        let (buf, write) = TempFileBuffer::new()?;
        let file = BufWriter::new(write);

        let (section_sender, section_receiver) = filebufferedchannel::channel(options.channel_size);
        let (sections_remote, sections_handle) =
            write_data(file, section_sender, frx).remote_handle();
        pool.spawn(sections_remote).expect("Couldn't spawn future.");
        (sections_handle, buf, section_receiver)
    };

    let processed_zooms: Vec<_> =
        std::iter::successors(Some(options.initial_zoom_size), |z| Some(z * 4))
            .take(options.max_zooms as usize)
            .map(|size| -> io::Result<_> {
                let (ftx, frx) = channel(options.channel_size);
                let (buf, write) = TempFileBuffer::new()?;
                let file = BufWriter::new(write);

                let (section_sender, section_receiver) =
                    filebufferedchannel::channel(options.channel_size);
                let (remote, handle) = write_data(file, section_sender, frx).remote_handle();
                pool.spawn(remote).expect("Couldn't spawn future.");
                let zoom_info = TempZoomInfo {
                    resolution: size,
                    data_write_future: Box::new(handle),
                    data: buf,
                    sections: section_receiver,
                };
                Ok((zoom_info, ftx))
            })
            .collect::<io::Result<_>>()?;
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
        },
    ))
}

#[cfg(test)]
mod tests {
    use byteordered::{ByteOrdered, Endianness};

    use super::*;
    use std::io::Cursor;

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
        let bufreader = BufReader::new(&mut cursor);
        let mut file = ByteOrdered::runtime(bufreader, Endianness::native());

        let magic = file.read_u32()?;
        if magic != CIR_TREE_MAGIC {
            panic!("Invalid file format: CIR_TREE_MAGIC does not match.");
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

        let mut blocks = vec![];
        crate::bbiread::search_overlapping_blocks(
            &mut file,
            Endianness::native(),
            0,
            0,
            MAX_BASES,
            &mut blocks,
        )?;

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
