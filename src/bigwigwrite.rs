use std::io::{self, Seek, SeekFrom, Write};
use std::io::{BufReader, BufWriter};
use std::fs::File;
use std::pin::Pin;
use std::vec::Vec;

use futures::try_join;
use futures::future::{Future, FutureExt, RemoteHandle};
use futures::channel::mpsc::{channel, Receiver};
use futures::executor::{block_on, ThreadPool};
use futures::sink::SinkExt;
use futures::stream::StreamExt;
use futures::task::SpawnExt;

use byteorder::{NativeEndian, WriteBytesExt};

use flate2::Compression;
use flate2::write::ZlibEncoder;

use serde::{Serialize, Deserialize};

use crate::chromvalues::ChromValues;
use crate::filebufferedchannel;
use crate::idmap::IdMap;
use crate::tell::Tell;
use crate::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};

use crate::bigwig::{Value, ZoomHeader, CHROM_TREE_MAGIC, CIR_TREE_MAGIC, BIGWIG_MAGIC_LTH};

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

type TempZoomInfo = (u32 /* resolution */, RemoteHandle<Result<usize, WriteGroupsError>> /* Temp file that contains data */, TempFileBuffer, filebufferedchannel::Receiver<Section> /* sections */);
type ZoomInfo = (u32 /* resolution */, File /* Temp file that contains data */, Box<dyn Iterator<Item=Section>> /* sections */);
const DEFAULT_ZOOM_SIZES: [u32; 11] = [10, 40, 160, 640, 2_560, 10_240, 40_960, 163_840, 655_360, 2_621_440, 10_485_760];

pub type ChromGroupRead = (
    Pin<Box<dyn Future<Output=Result<Summary, WriteGroupsError>> + Send>>,
    filebufferedchannel::Receiver<Section>,
    TempFileBuffer,
    Box<dyn Future<Output=Result<usize, WriteGroupsError>> + Send + Unpin>,
    Vec<TempZoomInfo>,
    (String, u32)
);

pub trait ChromGroupReadStreamingIterator {
    fn next(&mut self) -> Result<Option<ChromGroupRead>, WriteGroupsError>;
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

impl BigWigWrite {
    pub fn create_file(path: String) -> Self {
        BigWigWrite {
            path,
            options: BigWigWriteOptions {
                compress: true,
                items_per_slot: 1024,
                block_size: 256,
            }
        }
    }

    const MAX_ZOOM_LEVELS: usize = 10;

    pub fn write_groups<V>(&self, chrom_sizes: std::collections::HashMap<String, u32>, vals: V) -> Result<(), WriteGroupsError> where V : ChromGroupReadStreamingIterator + Send {
        let fp = File::create(self.path.clone())?;
        let mut file = BufWriter::new(fp);

        BigWigWrite::write_blank_headers(&mut file)?;

        let total_summary_offset = file.tell()?;
        file.write_all(&[0; 40])?;

        let full_data_offset = file.tell()?;

        // Total items
        // Unless we know the vals ahead of time, we can't estimate total sections ahead of time.
        // Even then simply doing "(vals.len() as u32 + ITEMS_PER_SLOT - 1) / ITEMS_PER_SLOT"
        // underestimates because sections are split by chrom too, not just size.
        // Skip for now, and come back when we write real header + summary.
        file.write_u32::<NativeEndian>(0)?;

        let pre_data = file.tell()?;
        // Write data to file and return
        let (chrom_ids, summary, mut file, raw_sections_iter, zoom_infos) = block_on(BigWigWrite::write_vals(vals, file))?;
        let data_size = file.tell()? - pre_data;
        let mut current_offset = pre_data;
        let sections_iter = raw_sections_iter.map(|mut section| {
            // TODO: this assumes that all the data is contiguous
            // This will fail if we ever space the sections in any way
            section.offset = current_offset;
            current_offset += section.size;
            section
        });
        //println!("Data size: {:?}", data_size);
        //println!("Sections: {:?}", total_sections);
        //println!("Summary: {:?}", summary);

        // Since the chrom tree is read before the index, we put this before the full data index
        // Therefore, there is a higher likelihood that the udc file will only need one read for chrom tree + full data index
        // Putting the chrom tree before the data also has a higher likelihood of being included with the beginning headers,
        // but requires us to know all the data ahead of time (when writing)
        let chrom_index_start = file.tell()?;
        BigWigWrite::write_chrom_tree(&mut file, chrom_sizes, &chrom_ids.get_map())?;

        let index_start = file.tell()?;
        let (nodes, levels, total_sections) = BigWigWrite::get_rtreeindex(sections_iter, &self.options);
        BigWigWrite::write_rtreeindex(&mut file, nodes, levels, total_sections, &self.options)?;

        let zoom_entries = BigWigWrite::write_zooms(&mut file, zoom_infos, data_size, &self.options)?;
        //println!("Zoom entries: {:?}", zoom_entries);
        let num_zooms = zoom_entries.len() as u16;
        //println!("Zooms: {:?}", num_zooms);

        // We *could* actually check the the real max size, but let's just assume at it's as large as the largest possible value
        // In most cases, I think this is the true max size (unless there is only one section and its less than ITEMS_PER_SLOT in size)
        let uncompress_buf_size = if self.options.compress {
            self.options.items_per_slot * (1 + 1 + 2 + 4 + 4 + 4 + 4 + 8 + 8)
        } else {
            0
        };

        file.seek(SeekFrom::Start(0))?;
        file.write_u32::<NativeEndian>(BIGWIG_MAGIC_LTH)?; // TODO: should really encode this with NativeEndian, since that is really what we do elsewhere
        file.write_u16::<NativeEndian>(4)?; // Actually 3, unsure what version 4 actually adds
        file.write_u16::<NativeEndian>(num_zooms)?;
        file.write_u64::<NativeEndian>(chrom_index_start)?;
        file.write_u64::<NativeEndian>(full_data_offset)?;
        file.write_u64::<NativeEndian>(index_start)?;
        file.write_u16::<NativeEndian>(0)?; // fieldCount
        file.write_u16::<NativeEndian>(0)?; // definedFieldCount
        file.write_u64::<NativeEndian>(0)?; // autoSQLOffset
        file.write_u64::<NativeEndian>(total_summary_offset)?;
        file.write_u32::<NativeEndian>(uncompress_buf_size)?;
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


        file.write_u32::<NativeEndian>(total_sections as u32)?;
        file.seek(SeekFrom::End(0))?;
        file.write_u32::<NativeEndian>(BIGWIG_MAGIC_LTH)?; // TODO: see above, should encode with NativeEndian

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

    async fn write_section(compress: bool, items_in_section: Vec<Value>, chrom_id: u32) -> io::Result<SectionData> {
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

    async fn process_group<I>(
        mut zooms_channels: Vec<futures::channel::mpsc::Sender<Pin<Box<dyn Future<Output=io::Result<SectionData>> + Send>>>>,
        mut ftx: futures::channel::mpsc::Sender<Pin<Box<dyn Future<Output=io::Result<SectionData>> + Send>>>,
        chrom_id: u32,
        options: BigWigWriteOptions,
        mut pool: ThreadPool,
        mut group: I,
        chrom: String,
        chrom_length: u32,
        ) -> Result<Summary, WriteGroupsError> 
        where I: ChromValues + Send {
        let num_zooms = DEFAULT_ZOOM_SIZES.len();
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
            // TODO: test these correctly fails
            if current_val.start > current_val.end {
                return Err(WriteGroupsError::InvalidInput(format!("Invalid bed graph: {} > {}", current_val.start, current_val.end)));
            }
            if current_val.start >= chrom_length {
                return Err(WriteGroupsError::InvalidInput(format!("Invalid bed graph: `{}` is greater than the chromosome ({}) length ({})", current_val.start, chrom, chrom_length)));
            }
            match group.peek() {
                None => (),
                Some(next_val) => {
                    if current_val.end > next_val.start {
                        return Err(WriteGroupsError::InvalidInput(format!(
                            "Invalid bed graph: overlapping values on chromosome {} at {}-{} and {}-{}",
                            chrom,
                            current_val.start,
                            current_val.end,
                            next_val.start,
                            next_val.end,
                        )));
                    }
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
                                chrom: chrom_id,
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
                                    chrom: chrom_id,
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
                let handle = pool.spawn_with_handle(BigWigWrite::write_section(options.compress, items, chrom_id)).expect("Couldn't spawn.");
                ftx.send(handle.boxed()).await.expect("Couldn't send");
            }
        }

        if !state_val.items.is_empty() {
            let handle = pool.spawn_with_handle(BigWigWrite::write_section(options.compress, state_val.items, chrom_id)).expect("Couldn't spawn.");
            ftx.send(handle.boxed()).await.expect("Couldn't send");
        }

        for (i, mut zoom_item) in state_val.zoom_items.into_iter().enumerate() {
            if let Some(zoom2) = zoom_item.live_info {
                debug_assert!(chrom_id == zoom2.chrom);
                zoom_item.records.push(ZoomRecord {
                    chrom: chrom_id,
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
    }

    async fn create_do_write<W: Write>(
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
                size: size,
            }).expect("Couldn't send section.");
            current_offset += size;
        }
        Ok(total)
    }

    pub fn read_group<I: 'static>(chrom: String, chrom_id: u32, chrom_length: u32, group: I, mut pool: ThreadPool, options: BigWigWriteOptions)
        -> io::Result<ChromGroupRead>
        where I: ChromValues + Send {
        let cloned_chrom = chrom.clone();

        let (ftx, frx) = channel::<_>(100);

        let (sections_future, buf, section_receiver) = {
            let (buf, write) = TempFileBuffer::new()?;
            let file = BufWriter::new(write);

            let (section_sender, section_receiver) = filebufferedchannel::channel::<Section>(200);
            let sections_future = BigWigWrite::create_do_write(file, section_sender, frx);
            (sections_future, buf, section_receiver)
        };

        let processed_zooms: Result<Vec<_>, _> = DEFAULT_ZOOM_SIZES.iter().map(|size| -> io::Result<_> {
            let (ftx, frx) = channel::<_>(100);
            let f = {
                let (buf, write) = TempFileBuffer::new()?;
                let file = BufWriter::new(write);

                let (section_sender, section_receiver) = filebufferedchannel::channel::<Section>(200);
                let file_future = BigWigWrite::create_do_write(file, section_sender, frx);

                (*size, file_future, buf, section_receiver)
            };
            Ok((f, ftx))
        }).collect();
        let (zooms_futures, zooms_channels): (Vec<_>, Vec<_>) = processed_zooms?.into_iter().unzip();

        let (sections_remote, sections_handle) = sections_future.remote_handle();
        let (zoom_infos, zoom_remotes): (Vec<_>, Vec<_>) = zooms_futures.into_iter().map(|(size, file_future, buf, section_iter)| {
            let (remote, handle) = file_future.remote_handle();
            ((size, handle, buf, section_iter), remote)
        }).unzip();

        for zoom_remote in zoom_remotes {
            pool.spawn(zoom_remote).expect("Couldn't spawn future.");
        }
        pool.spawn(sections_remote).expect("Couldn't spawn future.");
        let (f_remote, f_handle) = BigWigWrite::process_group(zooms_channels, ftx, chrom_id, options, pool.clone(), group, chrom, chrom_length).remote_handle();
        pool.spawn(f_remote).expect("Couldn't spawn future.");
        Ok((f_handle.boxed(), section_receiver, buf, Box::new(sections_handle), zoom_infos, (cloned_chrom, chrom_id)))    
    }

    async fn write_vals<V>(
        mut vals_iter: V,
        file: BufWriter<File>
    )
    -> Result<(
        IdMap<String>,
        Summary,
        BufWriter<File>,
        Box<dyn Iterator<Item=Section> + 'static>,
        Vec<ZoomInfo>
    ), WriteGroupsError>
    where V : ChromGroupReadStreamingIterator + Send {
        let mut section_iter: Vec<Box<dyn Iterator<Item=Section>>> = vec![];

        let mut zooms: Vec<(u32, Vec<Box<dyn Iterator<Item=Section>>>, TempFileBuffer, TempFileBufferWriter)> = DEFAULT_ZOOM_SIZES.iter().map(|size| -> io::Result<_> {
            let section_iter: Vec<Box<dyn Iterator<Item=Section>>> = vec![];

            let (buf, write) = TempFileBuffer::new()?;
            Ok((*size, section_iter, buf, write))
        }).collect::<Result<_, _>>()?;

        let mut raw_file = file.into_inner()?;

        let mut summary: Option<Summary> = None;

        let mut chrom_ids = IdMap::new();
        while let Some((summary_future, sections_receiver, mut sections_buf, sections_future, zoom_infos, (chrom, chrom_id))) = vals_iter.next()? {
            // We previously needed the chrom ids for the indices (read_group)
            // However, we still want the full map in order to write the full chrom tree
            // TODO: think about if this is the best way to do this (it's not) and alternatives
            // (maybe add a method to ChromGroupReadStreamingIterator `dissolve`)
            let real_id = chrom_ids.get_id(chrom);
            assert_eq!(real_id, chrom_id, "The chrom id passed from vals_iter does not match the expected chrom id.");

            // If we concurrently processing multiple chromosomes, the section buffer might have written some or all to a separate file
            // Switch that processing output to the real file
            sections_buf.switch(raw_file)?;
            // Await these futures to drive them to completion
            let (chrom_summary, _num_sections) = try_join!(summary_future, sections_future)?;
            section_iter.push(Box::new(sections_receiver.into_iter()));
            raw_file = sections_buf.await_file();

            for (i, (_size, future, buf, zoom_sections_receiver)) in zoom_infos.into_iter().enumerate() {
                let zoom = &mut zooms[i];
                assert_eq!(zoom.0, _size);
                let num_sections = future.await?;
                // TODO: do we need to `take` here?
                zoom.1.push(Box::new(zoom_sections_receiver.into_iter().take(num_sections)));
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
            let zoom_iter: Box<dyn Iterator<Item=Section> + 'static> = Box::new(zoom.1.into_iter().flat_map(|s| s));
            (zoom.0, zoom.2.await_file(), zoom_iter)
        }).collect();
        let section_iter: Box<dyn Iterator<Item=Section> + 'static> = Box::new(section_iter.into_iter().flat_map(|s| s));
        Ok((chrom_ids, summary_complete, BufWriter::new(raw_file), section_iter, zoom_infos))
    }

    fn write_zooms(mut file: &mut BufWriter<File>, zooms: Vec<ZoomInfo>, data_size: u64, options: &BigWigWriteOptions) -> io::Result<Vec<ZoomHeader>> {
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

        Ok(zoom_entries)
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


    const NODEHEADER_SIZE: u64 = 1 + 1 + 2;
    const NON_LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8;
    const LEAFNODE_SIZE: u64 = 4 + 4 + 4 + 4 + 8 + 8;

    fn calculate_offsets(mut index_offsets: &mut Vec<u64>, nodes: &RTreeChildren, level: usize) {
        match nodes {
            RTreeChildren::DataSections(_) => (),
            RTreeChildren::Nodes(children) => {
                index_offsets[level - 1] += BigWigWrite::NODEHEADER_SIZE;
                for child in children {
                    index_offsets[level - 1] += BigWigWrite::NON_LEAFNODE_SIZE;
                    BigWigWrite::calculate_offsets(&mut index_offsets, &child.children, level - 1);
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
        options: &BigWigWriteOptions
    ) -> io::Result<u64> {
        let non_leafnode_full_block_size: u64 = BigWigWrite::NODEHEADER_SIZE + BigWigWrite::NON_LEAFNODE_SIZE * options.block_size as u64;
        let leafnode_full_block_size: u64 = BigWigWrite::NODEHEADER_SIZE + BigWigWrite::LEAFNODE_SIZE * options.block_size as u64;
        debug_assert!(curr_level >= dest_level);
        let mut total_size = 0;
        if curr_level != dest_level {
            let mut next_offset_offset = 0;
            match nodes {
                RTreeChildren::DataSections(_) => panic!("Datasections found at level: {:?}", curr_level),
                RTreeChildren::Nodes(children) => {
                    for child in children {
                        next_offset_offset += BigWigWrite::write_tree(&mut file, &child.children, curr_level - 1, dest_level, childnode_offset + next_offset_offset, options)?;
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



    fn write_rtreeindex(
        file: &mut BufWriter<File>,
        nodes: RTreeChildren,
        levels: usize,
        section_count: u64,
        options: &BigWigWriteOptions
    ) -> io::Result<()> {
        let mut index_offsets: Vec<u64> = vec![0u64; levels as usize];

        BigWigWrite::calculate_offsets(&mut index_offsets, &nodes, levels);
        //println!("index Offsets: {:?}", index_offsets);

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
        //println!("Levels: {:?}", levels);
        //println!("Start of index: {}", next_offset);
        for level in (0..=levels).rev() {
            if level > 0 {
                next_offset += index_offsets[level - 1];
            }
            BigWigWrite::write_tree(file, &nodes, levels, level, next_offset, options)?;
            //println!("End of index level {}: {}", level, file.seek(SeekFrom::Current(0))?);
        }
        //println!("Total index size: {:?}", total_size);

        Ok(())
    }

}
