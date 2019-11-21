use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::pin::Pin;

use futures::executor::{block_on, ThreadPool};
use futures::future::{Future, FutureExt};
use futures::sink::SinkExt;
use futures::task::SpawnExt;

use byteorder::{NativeEndian, WriteBytesExt};

use crate::chromvalues::ChromValues;
use crate::tell::Tell;

use crate::bbiwrite::{
    encode_zoom_section, get_chromprocessing, get_rtreeindex, write_blank_headers,
    write_chrom_tree, write_rtreeindex, write_vals, write_zooms, BBIWriteOptions, ChromGroupRead,
    ChromGroupReadStreamingIterator, ChromProcessingInput, SectionData, WriteGroupsError,
};
use crate::bigwig::{BedEntry, Summary, Value, ZoomRecord, BIGBED_MAGIC};

pub struct BigBedWrite {
    pub path: String,
    pub options: BBIWriteOptions,
}

impl BigBedWrite {
    pub fn create_file(path: String) -> Self {
        BigBedWrite {
            path,
            options: BBIWriteOptions::default(),
        }
    }

    pub fn write_groups<V>(
        &self,
        chrom_sizes: HashMap<String, u32>,
        vals: V,
    ) -> Result<(), WriteGroupsError>
    where
        V: ChromGroupReadStreamingIterator + Send,
    {
        let fp = File::create(self.path.clone())?;
        let mut file = BufWriter::new(fp);

        write_blank_headers(&mut file)?;

        let total_summary_offset = file.tell()?;
        file.write_all(&[0; 40])?;

        // TODO: extra indices

        let full_data_offset = file.tell()?;

        // Total items
        // Unless we know the vals ahead of time, we can't estimate total sections ahead of time.
        // Even then simply doing "(vals.len() as u32 + ITEMS_PER_SLOT - 1) / ITEMS_PER_SLOT"
        // underestimates because sections are split by chrom too, not just size.
        // Skip for now, and come back when we write real header + summary.
        file.write_u64::<NativeEndian>(0)?;

        let pre_data = file.tell()?;
        // Write data to file and return
        let (chrom_ids, summary, mut file, raw_sections_iter, zoom_infos) =
            block_on(write_vals(vals, file, &self.options))?;
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
        let num_zooms = zoom_entries.len() as u16;

        // We *could* actually check the the real max size, but let's just assume at it's as large as the largest possible value
        // In most cases, I think this is the true max size (unless there is only one section and its less than ITEMS_PER_SLOT in size)
        let uncompress_buf_size = if self.options.compress {
            self.options.items_per_slot * (1 + 1 + 2 + 4 + 4 + 4 + 4 + 8 + 8)
        } else {
            0
        };

        file.seek(SeekFrom::Start(0))?;
        file.write_u32::<NativeEndian>(BIGBED_MAGIC)?; // TODO: should really encode this with NativeEndian, since that is really what we do elsewhere
        file.write_u16::<NativeEndian>(4)?; // Actually 3, unsure what version 4 actually adds
        file.write_u16::<NativeEndian>(num_zooms)?;
        file.write_u64::<NativeEndian>(chrom_index_start)?;
        file.write_u64::<NativeEndian>(full_data_offset)?;
        file.write_u64::<NativeEndian>(index_start)?;
        // TODO: actually write the correct value
        file.write_u16::<NativeEndian>(3)?; // fieldCount
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

        file.seek(SeekFrom::Start(full_data_offset))?;
        file.write_u64::<NativeEndian>(summary.total_items)?;
        file.seek(SeekFrom::End(0))?;
        file.write_u32::<NativeEndian>(BIGBED_MAGIC)?; // TODO: see above, should encode with NativeEndian

        Ok(())
    }

    async fn process_group<I>(
        mut zooms_channels: Vec<
            futures::channel::mpsc::Sender<
                Pin<Box<dyn Future<Output = io::Result<SectionData>> + Send>>,
            >,
        >,
        mut ftx: futures::channel::mpsc::Sender<
            Pin<Box<dyn Future<Output = io::Result<SectionData>> + Send>>,
        >,
        chrom_id: u32,
        options: BBIWriteOptions,
        pool: ThreadPool,
        mut group: I,
        chrom: String,
        chrom_length: u32,
    ) -> Result<Summary, WriteGroupsError>
    where
        I: ChromValues<BedEntry> + Send,
    {
        // While we do technically lose precision here by using the f32 in Value, we can reuse the same merge_into method
        struct ZoomItem {
            size: u32,
            live_info: Option<(ZoomRecord, u64)>,
            overlap: VecDeque<Value>,
            records: Vec<ZoomRecord>,
        }
        struct EntriesSection {
            items: Vec<BedEntry>,
            overlap: VecDeque<Value>,
            zoom_items: Vec<ZoomItem>,
        }

        let mut summary: Option<Summary> = None;

        let mut state_val = EntriesSection {
            items: Vec::with_capacity(options.items_per_slot as usize),
            overlap: VecDeque::new(),
            zoom_items: std::iter::successors(Some(options.initial_zoom_size), |z| Some(z * 4))
                .take(options.max_zooms as usize)
                .map(|size| ZoomItem {
                    size,
                    live_info: None,
                    overlap: VecDeque::new(),
                    records: Vec::with_capacity(options.items_per_slot as usize),
                })
                .collect(),
        };
        let mut total_items = 0;
        while let Some(current_val) = group.next()? {
            total_items += 1;
            // TODO: test these correctly fails
            if current_val.start > current_val.end {
                return Err(WriteGroupsError::InvalidInput(format!(
                    "Invalid bed: {} > {}",
                    current_val.start, current_val.end
                )));
            }
            if current_val.start >= chrom_length {
                return Err(WriteGroupsError::InvalidInput(format!(
                    "Invalid bed: `{}` is greater than the chromosome ({}) length ({})",
                    current_val.start, chrom, chrom_length
                )));
            }
            if let Some(next_val) = group.peek() {
                if current_val.start > next_val.start {
                    return Err(WriteGroupsError::InvalidInput(format!(
                        "Invalid bed: not sorted on chromosome {} at {}-{} (first) and {}-{} (second). Use sort -k1,1 -k2,2n to sort the bed before input.",
                        chrom,
                        current_val.start,
                        current_val.end,
                        next_val.start,
                        next_val.end,
                    )));
                }
            }

            let add_interval_to_summary =
                move |overlap: &mut VecDeque<Value>,
                      summary: &mut Option<Summary>,
                      item_start: u32,
                      item_end: u32,
                      next_start_opt: Option<u32>| {
                    // If any overlaps exists, it must be starting at the current start (else it would have to be after the current entry)
                    // If the overlap starts before, the entry wasn't correctly cut last iteration
                    debug_assert!(overlap
                        .front()
                        .map(|f| f.start == item_start)
                        .unwrap_or(true));

                    for o in overlap.iter_mut() {
                        o.value += 1.0;
                        if item_end < o.end {
                            o.end = item_end;
                            let end = o.end;
                            overlap.push_back(Value {
                                start: item_end,
                                end,
                                value: 1.0,
                            });
                            debug_assert_eq!(overlap.back().unwrap().start, item_end);
                            debug_assert_eq!(overlap.back().unwrap().end, end);
                            break;
                        }
                    }

                    debug_assert!(overlap.back().map(|o| o.end >= item_start).unwrap_or(true));

                    if overlap.back().map(|o| o.end).unwrap_or(item_start) == item_start {
                        overlap.push_back(Value {
                            start: item_start,
                            end: item_end,
                            value: 1.0,
                        })
                    }

                    let next_start = next_start_opt.unwrap_or(u32::max_value());

                    while overlap
                        .front()
                        .map(|f| f.start < next_start)
                        .unwrap_or(false)
                    {
                        let mut removed = overlap.pop_front().unwrap();
                        let (len, val) = if removed.end <= next_start {
                            (removed.end - removed.start, f64::from(removed.value))
                        } else {
                            let len = next_start - removed.start;
                            let val = f64::from(removed.value);
                            removed.start = next_start;
                            overlap.push_front(removed);
                            (len, val)
                        };

                        match summary {
                            None => {
                                *summary = Some(Summary {
                                    total_items: 0,
                                    bases_covered: u64::from(len),
                                    min_val: val,
                                    max_val: val,
                                    sum: f64::from(len) * val,
                                    sum_squares: f64::from(len) * val * val,
                                })
                            }
                            Some(summary) => {
                                summary.bases_covered += u64::from(len);
                                summary.min_val = summary.min_val.min(val);
                                summary.max_val = summary.max_val.max(val);
                                summary.sum += f64::from(len) * val;
                                summary.sum_squares += f64::from(len) * val * val;
                            }
                        }
                    }
                };

            add_interval_to_summary(
                &mut state_val.overlap,
                &mut summary,
                current_val.start,
                current_val.end,
                group.peek().map(|v| v.start),
            );

            for (i, zoom_item) in state_val.zoom_items.iter_mut().enumerate() {
                debug_assert_ne!(zoom_item.records.len(), options.items_per_slot as usize);

                let item_start = current_val.start;
                let item_end = current_val.end;
                let overlap = &mut zoom_item.overlap;
                for o in overlap.iter_mut() {
                    o.value += 1.0;
                    if item_end < o.end {
                        o.end = item_end;
                        let end = o.end;
                        overlap.push_back(Value {
                            start: item_end,
                            end,
                            value: 1.0,
                        });
                        debug_assert_eq!(overlap.back().unwrap().start, item_end);
                        debug_assert_eq!(overlap.back().unwrap().end, end);
                        break;
                    }
                }

                debug_assert!(overlap.back().map(|o| o.end >= item_start).unwrap_or(true));

                if overlap.back().map(|o| o.end).unwrap_or(item_start) == item_start {
                    overlap.push_back(Value {
                        start: item_start,
                        end: item_end,
                        value: 1.0,
                    })
                }

                let next_start = group.peek().map(|v| v.start).unwrap_or(u32::max_value());

                while overlap
                    .front()
                    .map(|f| f.start < next_start)
                    .unwrap_or(false)
                {
                    let mut removed = overlap.pop_front().unwrap();
                    let val = f64::from(removed.value);
                    let (removed_start, removed_end) = if removed.end <= next_start {
                        (removed.start, removed.end)
                    } else {
                        let start = removed.start;
                        removed.start = next_start;
                        overlap.push_front(removed);
                        (start, next_start)
                    };

                    let mut add_start = removed_start;
                    loop {
                        if add_start >= removed_end {
                            if group.peek().is_none() {
                                if let Some((mut zoom2, total_items)) = zoom_item.live_info.take() {
                                    zoom2.summary.total_items = total_items;
                                    zoom_item.records.push(zoom2);
                                }
                                if !zoom_item.records.is_empty() {
                                    let items = std::mem::replace(&mut zoom_item.records, vec![]);
                                    let handle = pool
                                        .spawn_with_handle(encode_zoom_section(
                                            options.compress,
                                            items,
                                        ))
                                        .expect("Couldn't spawn.");
                                    zooms_channels[i]
                                        .send(handle.boxed())
                                        .await
                                        .expect("Couln't send");
                                }
                            }
                            break;
                        }
                        let (zoom2, _) = zoom_item.live_info.get_or_insert((
                            ZoomRecord {
                                chrom: chrom_id,
                                start: add_start,
                                end: add_start,
                                summary: Summary {
                                    total_items: 0,
                                    bases_covered: 0,
                                    min_val: 1.0,
                                    max_val: 1.0,
                                    sum: 0.0,
                                    sum_squares: 0.0,
                                },
                            },
                            0,
                        ));
                        // The end of zoom record
                        let next_end = zoom2.start + zoom_item.size;
                        // End of bases that we could add
                        let add_end = std::cmp::min(next_end, removed_end);
                        // If the last zoom ends before this value starts, we don't add anything
                        if add_end >= add_start {
                            let added_bases = add_end - add_start;
                            zoom2.end = add_end;
                            zoom2.summary.total_items += 1; // XXX
                            zoom2.summary.bases_covered += u64::from(added_bases);
                            zoom2.summary.min_val = zoom2.summary.min_val.min(val);
                            zoom2.summary.max_val = zoom2.summary.max_val.max(val);
                            zoom2.summary.sum += f64::from(added_bases) * val;
                            zoom2.summary.sum_squares += f64::from(added_bases) * val * val;
                        }
                        // If we made it to the end of the zoom (whether it was because the zoom ended before this value started,
                        // or we added to the end of the zoom), then write this zooms to the current section
                        if add_end == next_end {
                            zoom_item.records.push(
                                zoom_item
                                    .live_info
                                    .take()
                                    .map(|(mut zoom_item, total_items)| {
                                        zoom_item.summary.total_items = total_items;
                                        zoom_item
                                    })
                                    .unwrap(),
                            );
                        }
                        // Set where we would start for next time
                        add_start = std::cmp::max(add_end, removed_start);
                        // Write section if full
                        if zoom_item.records.len() == options.items_per_slot as usize {
                            let items = std::mem::replace(&mut zoom_item.records, vec![]);
                            let handle = pool
                                .spawn_with_handle(encode_zoom_section(options.compress, items))
                                .expect("Couldn't spawn.");
                            zooms_channels[i]
                                .send(handle.boxed())
                                .await
                                .expect("Couln't send");
                        }
                    }
                }

                debug_assert_ne!(zoom_item.records.len(), options.items_per_slot as usize);
            }
            state_val.items.push(current_val);
            if group.peek().is_none() || state_val.items.len() >= options.items_per_slot as usize {
                let items = std::mem::replace(&mut state_val.items, vec![]);
                let handle = pool
                    .spawn_with_handle(encode_section(options.compress, items, chrom_id))
                    .expect("Couldn't spawn.");
                ftx.send(handle.boxed()).await.expect("Couldn't send");
            }
        }

        debug_assert!(state_val.items.is_empty());
        for zoom_item in state_val.zoom_items.iter_mut() {
            debug_assert!(zoom_item.live_info.is_none());
            debug_assert!(zoom_item.records.is_empty());
        }

        let mut summary_complete = match summary {
            None => Summary {
                total_items: 0,
                bases_covered: 0,
                min_val: 0.0,
                max_val: 0.0,
                sum: 0.0,
                sum_squares: 0.0,
            },
            Some(summary) => summary,
        };
        summary_complete.total_items = total_items;
        Ok(summary_complete)
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
    /// The futures that are returned are only handles to remote futures that are spawned immediately on `pool`.
    pub fn begin_processing_chrom<I: 'static>(
        chrom: String,
        chrom_id: u32,
        chrom_length: u32,
        group: I,
        mut pool: ThreadPool,
        options: BBIWriteOptions,
    ) -> io::Result<ChromGroupRead>
    where
        I: ChromValues<BedEntry> + Send,
    {
        let (
            ChromProcessingInput {
                zooms_channels,
                ftx,
            },
            processing_output,
        ) = get_chromprocessing(&mut pool, &options)?;

        let (f_remote, f_handle) = BigBedWrite::process_group(
            zooms_channels,
            ftx,
            chrom_id,
            options,
            pool.clone(),
            group,
            chrom,
            chrom_length,
        )
        .remote_handle();
        pool.spawn(f_remote).expect("Couldn't spawn future.");
        let read = ChromGroupRead {
            summary_future: f_handle.boxed(),
            processing_output,
        };
        Ok(read)
    }
}

async fn encode_section(
    compress: bool,
    items_in_section: Vec<BedEntry>,
    chrom_id: u32,
) -> io::Result<SectionData> {
    use libdeflater::{CompressionLvl, Compressor};

    let mut bytes: Vec<u8> = vec![];

    let start = items_in_section[0].start;
    let end = items_in_section[items_in_section.len() - 1].end;

    for item in items_in_section.iter() {
        bytes.write_u32::<NativeEndian>(chrom_id)?;
        bytes.write_u32::<NativeEndian>(item.start)?;
        bytes.write_u32::<NativeEndian>(item.end)?;
        bytes.write_all(item.rest.as_bytes())?;
        bytes.write_all(&[b'\0'])?;
    }

    let out_bytes = if compress {
        let mut compressor = Compressor::new(CompressionLvl::default());
        let max_sz = compressor.zlib_compress_bound(bytes.len());
        let mut compressed_data = vec![0; max_sz];
        let actual_sz = compressor
            .zlib_compress(&bytes, &mut compressed_data)
            .unwrap();
        compressed_data.resize(actual_sz, 0);
        compressed_data
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
