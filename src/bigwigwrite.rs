use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::fs::File;
use std::pin::Pin;
use std::vec::Vec;

use futures::try_join;
use futures::future::{Either, Future, FutureExt, RemoteHandle};
use futures::channel::mpsc::{channel, Receiver};
use futures::executor::{block_on, ThreadPool};
use futures::sink::SinkExt;
use futures::stream::StreamExt;
use futures::task::SpawnExt;

use byteorder::{NativeEndian, WriteBytesExt};

use crate::chromvalues::ChromValues;
use crate::filebufferedchannel;
use crate::idmap::IdMap;
use crate::tell::Tell;
use crate::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};

use crate::bigwig::{Value, BIGWIG_MAGIC_LTH};
use crate::bbiwrite::{
    write_blank_headers,
    encode_section,
    encode_zoom_section,
    write_chrom_tree,
    Section,
    SectionData,
    ZoomRecord,
    get_rtreeindex,
    write_rtreeindex,
    BBIWriteOptions,
    ZoomInfo,
    write_zooms
};

#[derive(Debug)]
pub struct Summary {
    bases_covered: u64,
    min_val: f64,
    max_val: f64,
    sum: f64,
    sum_squares: f64,
}

type TempZoomInfo = (u32 /* resolution */, RemoteHandle<Result<usize, WriteGroupsError>> /* Temp file that contains data */, TempFileBuffer<TempFileBufferWriter<File>>, filebufferedchannel::Receiver<Section> /* sections */);
const DEFAULT_ZOOM_SIZES: [u32; 11] = [10, 40, 160, 640, 2_560, 10_240, 40_960, 163_840, 655_360, 2_621_440, 10_485_760];

pub type ChromGroupRead = (
    Pin<Box<dyn Future<Output=Result<Summary, WriteGroupsError>> + Send>>,
    filebufferedchannel::Receiver<Section>,
    TempFileBuffer<File>,
    Box<dyn Future<Output=Result<usize, WriteGroupsError>> + Send + Unpin>,
    Vec<TempZoomInfo>,
);

pub trait ChromGroupReadStreamingIterator {
    fn next(&mut self) -> Result<Option<Either<ChromGroupRead, IdMap<String>>>, WriteGroupsError>;
}

pub struct BigWigWrite {
    pub path: String,
    pub options: BBIWriteOptions,
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
            options: BBIWriteOptions {
                compress: true,
                items_per_slot: 1024,
                block_size: 256,
            }
        }
    }

    pub fn write_groups<V>(&self, chrom_sizes: std::collections::HashMap<String, u32>, vals: V) -> Result<(), WriteGroupsError> where V : ChromGroupReadStreamingIterator + Send {
        let fp = File::create(self.path.clone())?;
        let mut file = BufWriter::new(fp);

        write_blank_headers(&mut file)?;

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

        // Since the chrom tree is read before the index, we put this before the full data index
        // Therefore, there is a higher likelihood that the udc file will only need one read for chrom tree + full data index
        // Putting the chrom tree before the data also has a higher likelihood of being included with the beginning headers,
        // but requires us to know all the data ahead of time (when writing)
        let chrom_index_start = file.tell()?;
        write_chrom_tree(&mut file, chrom_sizes, &chrom_ids.get_map())?;

        let index_start = file.tell()?;
        let (nodes, levels, total_sections) = get_rtreeindex(sections_iter, &self.options);
        write_rtreeindex(&mut file, nodes, levels, total_sections, &self.options)?;

        let zoom_entries = write_zooms(&mut file, zoom_infos, data_size, &self.options)?;
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

    async fn process_group<I>(
        mut zooms_channels: Vec<futures::channel::mpsc::Sender<Pin<Box<dyn Future<Output=io::Result<SectionData>> + Send>>>>,
        mut ftx: futures::channel::mpsc::Sender<Pin<Box<dyn Future<Output=io::Result<SectionData>> + Send>>>,
        chrom_id: u32,
        options: BBIWriteOptions,
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
                                let handle = pool.spawn_with_handle(encode_zoom_section(options.compress, items)).expect("Couldn't spawn.");
                                zooms_channels[i].send(handle.boxed()).await.expect("Couln't send");
                            }
                        }
                    }
                }
            }
            state_val.items.push(current_val);
            if state_val.items.len() >= options.items_per_slot as usize {
                let items = std::mem::replace(&mut state_val.items, vec![]);
                let handle = pool.spawn_with_handle(encode_section(options.compress, items, chrom_id)).expect("Couldn't spawn.");
                ftx.send(handle.boxed()).await.expect("Couldn't send");
            }
        }

        if !state_val.items.is_empty() {
            let handle = pool.spawn_with_handle(encode_section(options.compress, state_val.items, chrom_id)).expect("Couldn't spawn.");
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
                let handle = pool.spawn_with_handle(encode_zoom_section(options.compress, items)).expect("Couldn't spawn.");
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

    /// This converts a ChromValues (streaming iterator) to a ChromGroupRead.
    /// This is a separate function so this can techincally be run for mulitple chromosomes simulatenously.
    /// This is heavily multi-threaded using Futures. A brief summary:
    /// - All reading from the ChromValues is done in a single future (process_group). This futures in charge of keeping track of sections (and zoom data).
    ///   When a section is full, a Future is created to byte-encode the data and compress it (if compression is on). The same is true for zoom sections.
    ///   This is the most CPU-intensive part of the entire write.
    /// - The section futures are sent (in order) by channel to a separate future for the sole purpose of writing the (maybe compressed) section data to a `TempFileBuffer`.
    ///   The data is written to a temporary file, since this may be happening in parallel (where we don't know the real file offset of these sections).
    ///   `TempFileBuffer` allows "switching" to the real file once it's available (on the read side), so if the real file is available, I/O ops are not duplicated.
    ///   Once the section data is written to the file, the file offset data is stored in a `FileBufferedChannel`. This is needed for the index.
    ///   All of this is done for zoom sections too.
    ///
    /// The returned ChromGroupRead contains:
    /// - The first future. This is important to be awaited for summary data and any errors.
    /// - The `FileBufferedChannel` with section offset data, for the index.
    /// - The `TempFileBuffer` to get the written section data
    /// - The second future. This is important to be awaited for any errors and the total number of sections wrriten.
    /// - A Vec of essentially the last three items, but for each zoom.
    /// 
    /// The futures that are returned are only handles to remote futures that are spawned immediately on `pool`.
    pub fn begin_processing_chrom<I: 'static>(chrom: String, chrom_id: u32, chrom_length: u32, group: I, mut pool: ThreadPool, options: BBIWriteOptions)
        -> io::Result<ChromGroupRead>
        where I: ChromValues + Send {

        let (ftx, frx) = channel::<_>(100);

        let (sections_handle, buf, section_receiver) = {
            let (buf, write) = TempFileBuffer::new()?;
            let file = BufWriter::new(write);

            let (section_sender, section_receiver) = filebufferedchannel::channel::<Section>(200);
            let (sections_remote, sections_handle) = BigWigWrite::create_do_write(file, section_sender, frx).remote_handle();
            pool.spawn(sections_remote).expect("Couldn't spawn future.");
            (sections_handle, buf, section_receiver)
        };

        let processed_zooms: Vec<_> = DEFAULT_ZOOM_SIZES.iter().map(|size| -> io::Result<_> {
            let (ftx, frx) = channel::<_>(100);
            let (buf, write) = TempFileBuffer::new()?;
            let file = BufWriter::new(write);

            let (section_sender, section_receiver) = filebufferedchannel::channel::<Section>(200);
            let (remote, handle) = BigWigWrite::create_do_write(file, section_sender, frx).remote_handle();
            pool.spawn(remote).expect("Couldn't spawn future.");
            Ok(((*size, handle, buf, section_receiver), ftx))
        }).collect::<io::Result<_>>()?;
        let (zoom_infos, zooms_channels): (Vec<_>, Vec<_>) = processed_zooms.into_iter().unzip();

        let (f_remote, f_handle) = BigWigWrite::process_group(zooms_channels, ftx, chrom_id, options, pool.clone(), group, chrom, chrom_length).remote_handle();
        pool.spawn(f_remote).expect("Couldn't spawn future.");
        Ok((f_handle.boxed(), section_receiver, buf, Box::new(sections_handle), zoom_infos))
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

        // Zooms have to be double-buffered: first because chroms could be processed in parallel and second because we don't know the offset of each zoom immediately
        let mut zooms: Vec<(u32, Vec<Box<dyn Iterator<Item=Section>>>, TempFileBuffer<File>, Option<TempFileBufferWriter<File>>)> = DEFAULT_ZOOM_SIZES.iter().map(|size| -> io::Result<_> {
            let section_iter: Vec<Box<dyn Iterator<Item=Section>>> = vec![];

            let (buf, write): (TempFileBuffer<File>, TempFileBufferWriter<File>) = TempFileBuffer::new()?;
            Ok((*size, section_iter, buf, Some(write)))
        }).collect::<Result<Vec<(u32, Vec<Box<dyn Iterator<Item=Section>>>, TempFileBuffer<File>, Option<TempFileBufferWriter<File>>)>, _>>()?;

        let mut raw_file = file.into_inner()?;

        let mut summary: Option<Summary> = None;

        let chrom_ids = loop {
            let next = vals_iter.next()?;
            match next {
                Some(Either::Left((summary_future, sections_receiver, mut sections_buf, sections_future, mut zoom_infos))) => {
                    // If we concurrently processing multiple chromosomes, the section buffer might have written some or all to a separate file
                    // Switch that processing output to the real file
                    sections_buf.switch(raw_file);
                    for (i, (size, _, buf, _)) in zoom_infos.iter_mut().enumerate() {
                        let zoom = &mut zooms[i];
                        let writer = zoom.3.take().unwrap();
                        assert_eq!(zoom.0, *size);
                        buf.switch(writer);
                    }
                    
                    // Await these futures to drive them to completion
                    let (chrom_summary, _num_sections) = try_join!(summary_future, sections_future)?;
                    section_iter.push(Box::new(sections_receiver.into_iter()));
                    raw_file = sections_buf.await_real_file();

                    for (i, (_, future, buf, zoom_sections_receiver)) in zoom_infos.into_iter().enumerate() {
                        let zoom = &mut zooms[i];
                        let _num_sections = future.await?;
                        zoom.1.push(Box::new(zoom_sections_receiver.into_iter()));
                        zoom.3.replace(buf.await_real_file());
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
                },
                Some(Either::Right(chrom_ids)) => break chrom_ids,
                None => unreachable!(),
            }
        };

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
            let closed_file = zoom.2.await_temp_file();
            (zoom.0, closed_file, zoom_iter)
        }).collect();
        let section_iter: Box<dyn Iterator<Item=Section> + 'static> = Box::new(section_iter.into_iter().flat_map(|s| s));
        Ok((chrom_ids, summary_complete, BufWriter::new(raw_file), section_iter, zoom_infos))
    }

}
