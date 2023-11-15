use std::collections::HashMap;
use std::ffi::CString;
use std::fs::File;
use std::io::{self, BufWriter, Write};

use futures::future::FutureExt;
use futures::sink::SinkExt;

use byteorder::{NativeEndian, WriteBytesExt};
use tokio::runtime::{Handle, Runtime};
use tokio::task;

use crate::utils::chromvalues::ChromValues;
use crate::utils::indexlist::IndexList;
use crate::utils::tell::Tell;
use crate::{write_info, ChromData, ChromProcessingInputSectionChannel};

use crate::bbi::{BedEntry, Summary, Value, ZoomRecord, BIGBED_MAGIC};
use crate::bbiwrite::{
    self, encode_zoom_section, get_rtreeindex, write_blank_headers, write_chrom_tree,
    write_rtreeindex, write_zooms, BBIWriteOptions, ProcessChromError, SectionData,
};

/// The struct used to write a bigBed file
pub struct BigBedWrite {
    pub path: String,
    pub options: BBIWriteOptions,
    pub autosql: Option<String>,
}

impl BigBedWrite {
    pub fn create_file(path: String) -> Self {
        BigBedWrite {
            path,
            options: BBIWriteOptions::default(),
            autosql: None,
        }
    }

    pub fn write<
        Values: ChromValues<Value = BedEntry> + Send + 'static,
        V: ChromData<Values = Values>,
    >(
        self,
        chrom_sizes: HashMap<String, u32>,
        vals: V,
        runtime: Runtime,
    ) -> Result<(), ProcessChromError<Values::Error>> {
        let runtime_handle = runtime.handle();

        let fp = File::create(self.path.clone())?;
        let mut file = BufWriter::new(fp);

        write_blank_headers(&mut file)?;

        let autosql_offset = file.tell()?;
        let autosql = self
            .autosql
            .clone()
            .unwrap_or_else(|| crate::bed::autosql::BED3.to_string());
        let autosql = CString::new(autosql.into_bytes()).map_err(|_| {
            ProcessChromError::InvalidInput("Invalid autosql: null byte in string".to_owned())
        })?;
        file.write_all(autosql.as_bytes_with_nul())?;

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

        let local = task::LocalSet::new();
        let output = runtime.block_on(local.run_until(bbiwrite::write_vals(
            vals,
            file,
            self.options,
            BigBedWrite::process_chrom,
            runtime_handle.clone(),
            chrom_sizes.clone(),
        )));
        let (chrom_ids, summary, mut file, raw_sections_iter, zoom_infos, uncompress_buf_size) =
            output?;

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
        write_chrom_tree(&mut file, chrom_sizes, &chrom_ids.get_map())?;

        let index_start = file.tell()?;
        let (nodes, levels, total_sections) = get_rtreeindex(sections_iter, self.options);
        write_rtreeindex(&mut file, nodes, levels, total_sections, self.options)?;

        let zoom_entries = write_zooms(&mut file, zoom_infos, data_size, self.options)?;
        let num_zooms = zoom_entries.len() as u16;

        write_info(
            &mut file,
            BIGBED_MAGIC,
            num_zooms,
            chrom_index_start,
            full_data_offset,
            index_start,
            // TODO: actually write the correct values for the following
            3,
            3,
            autosql_offset,
            total_summary_offset,
            uncompress_buf_size,
            zoom_entries,
            summary,
            summary.total_items,
        )?;

        Ok(())
    }

    async fn process_chrom<I>(
        zooms_channels: Vec<(u32, ChromProcessingInputSectionChannel)>,
        mut ftx: ChromProcessingInputSectionChannel,
        chrom_id: u32,
        options: BBIWriteOptions,
        runtime: Handle,
        mut group: I,
        chrom: String,
        chrom_length: u32,
    ) -> Result<Summary, ProcessChromError<I::Error>>
    where
        I: ChromValues<Value = BedEntry> + Send,
    {
        // While we do technically lose precision here by using the f32 in Value, we can reuse the same merge_into method
        struct ZoomItem {
            size: u32,
            live_info: Option<(ZoomRecord, u64)>,
            overlap: IndexList<Value>,
            records: Vec<ZoomRecord>,
            channel: ChromProcessingInputSectionChannel,
        }
        struct EntriesSection {
            items: Vec<BedEntry>,
            overlap: IndexList<Value>,
            zoom_items: Vec<ZoomItem>,
        }

        let mut summary: Option<Summary> = None;

        let mut state_val = EntriesSection {
            items: Vec::with_capacity(options.items_per_slot as usize),
            overlap: IndexList::new(),
            zoom_items: zooms_channels
                .into_iter()
                .map(|(size, channel)| ZoomItem {
                    size,
                    live_info: None,
                    overlap: IndexList::new(),
                    records: Vec::with_capacity(options.items_per_slot as usize),
                    channel,
                })
                .collect(),
        };
        let mut total_items = 0;
        while let Some(current_val) = group.next() {
            let current_val = match current_val {
                Ok(v) => v,
                Err(e) => return Err(ProcessChromError::SourceError(e)),
            };
            total_items += 1;
            // TODO: test these correctly fails
            if current_val.start > current_val.end {
                return Err(ProcessChromError::InvalidInput(format!(
                    "Invalid bed: {} > {}",
                    current_val.start, current_val.end
                )));
            }
            if current_val.start >= chrom_length {
                return Err(ProcessChromError::InvalidInput(format!(
                    "Invalid bed: `{}` is greater than the chromosome ({}) length ({})",
                    current_val.start, chrom, chrom_length
                )));
            }
            if let Some(Ok(next_val)) = group.peek() {
                if current_val.start > next_val.start {
                    return Err(ProcessChromError::InvalidInput(format!(
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
                move |overlap: &mut IndexList<Value>,
                      summary: &mut Option<Summary>,
                      item_start: u32,
                      item_end: u32,
                      next_start_opt: Option<u32>| {
                    // If any overlaps exists, it must be starting at the current start (else it would have to be after the current entry)
                    // If the overlap starts before, the entry wasn't correctly cut last iteration
                    debug_assert!(overlap
                        .head()
                        .map(|f| f.start == item_start)
                        .unwrap_or(true));

                    // For each item in `overlap` that overlaps the current
                    // item, add `1` to the value.
                    let mut index = overlap.head_index();
                    while let Some(i) = index {
                        match overlap.get_mut(i) {
                            None => break,
                            Some(o) => {
                                o.value += 1.0;
                                if item_end < o.end {
                                    let value = o.value - 1.0;
                                    let end = o.end;
                                    o.end = item_end;
                                    overlap.insert_after(
                                        i,
                                        Value {
                                            start: item_end,
                                            end,
                                            value,
                                        },
                                    );
                                    break;
                                }
                                index = overlap.next_index(i);
                            }
                        }
                    }

                    debug_assert!(overlap.tail().map(|o| o.end >= item_start).unwrap_or(true));

                    if overlap.tail().map(|o| o.end).unwrap_or(item_start) == item_start {
                        overlap.push_back(Value {
                            start: item_start,
                            end: item_end,
                            value: 1.0,
                        });
                    }

                    let next_start = next_start_opt.unwrap_or(u32::max_value());

                    while overlap
                        .head()
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
                group.peek().and_then(|v| v.ok()).map(|v| v.start),
            );

            for zoom_item in state_val.zoom_items.iter_mut() {
                debug_assert_ne!(zoom_item.records.len(), options.items_per_slot as usize);

                let item_start = current_val.start;
                let item_end = current_val.end;
                let overlap = &mut zoom_item.overlap;

                // For each item in `overlap` that overlaps the current
                // item, add `1` to the value.
                let mut index = overlap.head_index();
                while let Some(i) = index {
                    match overlap.get_mut(i) {
                        None => break,
                        Some(o) => {
                            o.value += 1.0;
                            if item_end < o.end {
                                let value = o.value - 1.0;
                                let end = o.end;
                                o.end = item_end;
                                overlap.insert_after(
                                    i,
                                    Value {
                                        start: item_end,
                                        end,
                                        value,
                                    },
                                );
                                break;
                            }
                            index = overlap.next_index(i);
                        }
                    }
                }

                debug_assert!(overlap.tail().map(|o| o.end >= item_start).unwrap_or(true));

                if overlap.tail().map(|o| o.end).unwrap_or(item_start) == item_start {
                    overlap.push_back(Value {
                        start: item_start,
                        end: item_end,
                        value: 1.0,
                    });
                }

                let next_start = group
                    .peek()
                    .and_then(|v| v.ok())
                    .map(|v| v.start)
                    .unwrap_or(u32::max_value());

                while overlap
                    .head()
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
                                    let items = std::mem::take(&mut zoom_item.records);
                                    let handle = runtime
                                        .spawn(encode_zoom_section(options.compress, items))
                                        .map(|f| f.unwrap());
                                    zoom_item
                                        .channel
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
                            let items = std::mem::take(&mut zoom_item.records);
                            let handle = runtime
                                .spawn(encode_zoom_section(options.compress, items))
                                .map(|f| f.unwrap());
                            zoom_item
                                .channel
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
                let items = std::mem::replace(
                    &mut state_val.items,
                    Vec::with_capacity(options.items_per_slot as usize),
                );
                let handle = runtime
                    .spawn(encode_section(options.compress, items, chrom_id))
                    .map(|f| f.unwrap());
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
}

async fn encode_section(
    compress: bool,
    items_in_section: Vec<BedEntry>,
    chrom_id: u32,
) -> io::Result<(SectionData, usize)> {
    use libdeflater::{CompressionLvl, Compressor};

    let mut bytes = Vec::with_capacity(items_in_section.len() * 30);

    let start = items_in_section[0].start;
    let end = items_in_section[items_in_section.len() - 1].end;

    // FIXME: Each of these calls end up calling `Vec::reserve`
    // We could instead use a `Cursor<&mut [u8]>`, but we would need to be a bit
    // more careful here around safety
    for item in items_in_section.iter() {
        bytes.write_u32::<NativeEndian>(chrom_id)?;
        bytes.write_u32::<NativeEndian>(item.start)?;
        bytes.write_u32::<NativeEndian>(item.end)?;
        bytes.write_all(item.rest.as_bytes())?;
        bytes.write_all(&[b'\0'])?;
    }

    let (out_bytes, uncompress_buf_size) = if compress {
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
            chrom: chrom_id,
            start,
            end,
            data: out_bytes,
        },
        uncompress_buf_size,
    ))
}
