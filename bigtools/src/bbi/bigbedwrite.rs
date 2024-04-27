use std::collections::HashMap;
use std::error::Error;
use std::ffi::CString;
use std::fs::File;
use std::io::{self, BufWriter, Write};

use futures::future::FutureExt;
use futures::sink::SinkExt;
use futures::Future;

use byteorder::{NativeEndian, WriteBytesExt};
use tokio::runtime::{Handle, Runtime};

use crate::utils::chromvalues::ChromValues;
use crate::utils::indexlist::IndexList;
use crate::utils::tell::Tell;
use crate::{write_info, ChromData, ChromProcessingInputSectionChannel};

use crate::bbi::{BedEntry, Summary, Value, ZoomRecord, BIGBED_MAGIC};
use crate::bbiwrite::{
    self, encode_zoom_section, write_blank_headers, write_zooms, BBIWriteOptions,
    ProcessChromError, SectionData,
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

    fn write_pre<E: Error>(
        file: &mut BufWriter<File>,
        autosql: &Option<String>,
    ) -> Result<(u64, u64, u64, u64), ProcessChromError<E>> {
        write_blank_headers(file)?;

        let autosql_offset = file.tell()?;
        let autosql = autosql
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

        Ok((
            autosql_offset,
            total_summary_offset,
            full_data_offset,
            pre_data,
        ))
    }

    /// Write the values from `V` as a bigWig. Will utilize the provided runtime for encoding values and for reading through the values (potentially parallelized by chromosome).
    pub fn write<
        Values: ChromValues<Value = BedEntry> + Send + 'static,
        V: ChromData<Values = Values>,
    >(
        self,
        chrom_sizes: HashMap<String, u32>,
        vals: V,
        runtime: Runtime,
    ) -> Result<(), ProcessChromError<Values::Error>> {
        let process_chrom = |zooms_channels: Vec<(u32, ChromProcessingInputSectionChannel)>,
                             ftx: ChromProcessingInputSectionChannel,
                             chrom_id: u32,
                             options: BBIWriteOptions,
                             runtime: Handle,
                             chrom_values: Values,
                             chrom: String,
                             chrom_length: u32| {
            let fut = BigBedWrite::process_chrom(
                zooms_channels,
                ftx,
                chrom_id,
                options,
                runtime.clone(),
                chrom_values,
                chrom,
                chrom_length,
            );
            runtime.spawn(fut).map(|f| f.unwrap())
        };
        self.write_internal(chrom_sizes, vals, runtime, process_chrom)
    }

    /// Write the values from `V` as a bigWig. Will utilize the provided runtime for encoding values, but will read through values on the current thread.
    pub fn write_singlethreaded<
        Values: ChromValues<Value = BedEntry>,
        V: ChromData<Values = Values>,
    >(
        self,
        chrom_sizes: HashMap<String, u32>,
        vals: V,
        runtime: Runtime,
    ) -> Result<(), ProcessChromError<Values::Error>> {
        self.write_internal(chrom_sizes, vals, runtime, BigBedWrite::process_chrom)
    }

    fn write_internal<
        Values: ChromValues<Value = BedEntry>,
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
        self,
        chrom_sizes: HashMap<String, u32>,
        vals: V,
        runtime: Runtime,
        process_chrom: G,
    ) -> Result<(), ProcessChromError<Values::Error>> {
        let fp = File::create(self.path.clone())?;
        let mut file = BufWriter::new(fp);

        let (autosql_offset, total_summary_offset, full_data_offset, pre_data) =
            BigBedWrite::write_pre(&mut file, &self.autosql)?;

        let output = bbiwrite::write_vals(
            vals,
            file,
            self.options,
            process_chrom,
            runtime,
            chrom_sizes.clone(),
        );
        let (chrom_ids, summary, mut file, raw_sections_iter, zoom_infos, uncompress_buf_size) =
            output?;

        let chrom_ids = chrom_ids.get_map();
        let (data_size, chrom_index_start, index_start, _total_sections) = bbiwrite::write_mid(
            &mut file,
            pre_data,
            raw_sections_iter,
            chrom_sizes,
            &chrom_ids,
            self.options,
        )?;

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
            // In bigWigs, this is total sections, but total items in bigBeds
            summary.total_items,
        )?;

        Ok(())
    }

    /// Write the values from `V` as a bigBed. Will utilize the provided runtime for encoding values and for reading through the values (potentially parallelized by chromosome).
    /// This will take two passes on the provided values: first to write the values themselves, then the zooms. This is beneficial over `write` on smaller files, where the encoding of
    /// high resolution zooms takes up a substantial portion of total processing time.
    pub fn write_multipass<
        Values: ChromValues<Value = BedEntry> + Send + 'static,
        V: ChromData<Values = Values>,
    >(
        self,
        make_vals: impl Fn() -> Result<V, ProcessChromError<Values::Error>>,
        chrom_sizes: HashMap<String, u32>,
        runtime: Runtime,
    ) -> Result<(), ProcessChromError<Values::Error>> {
        let fp = File::create(self.path.clone())?;
        let mut file = BufWriter::new(fp);

        let (autosql_offset, total_summary_offset, full_data_offset, pre_data) =
            BigBedWrite::write_pre(&mut file, &self.autosql)?;

        let vals = make_vals()?;

        let runtime_handle = runtime.handle();

        let process_chrom = |ftx: ChromProcessingInputSectionChannel,
                             chrom_id: u32,
                             options: BBIWriteOptions,
                             runtime: Handle,
                             chrom_values: Values,
                             chrom: String,
                             chrom_length: u32| {
            let fut = BigBedWrite::process_chrom_no_zooms(
                ftx,
                chrom_id,
                options,
                runtime,
                chrom_values,
                chrom,
                chrom_length,
            );
            let (fut, handle) = fut.remote_handle();
            runtime_handle.spawn(fut);
            handle
        };

        let output = bbiwrite::write_vals_no_zoom(
            vals,
            file,
            self.options,
            process_chrom,
            &runtime,
            chrom_sizes.clone(),
        );
        let (chrom_ids, summary, zoom_counts, mut file, raw_sections_iter, mut uncompress_buf_size) =
            output?;

        let chrom_ids = chrom_ids.get_map();
        let (data_size, chrom_index_start, index_start, _total_sections) = bbiwrite::write_mid(
            &mut file,
            pre_data,
            raw_sections_iter,
            chrom_sizes,
            &chrom_ids,
            self.options,
        )?;

        let vals = make_vals()?;

        let output = bbiwrite::write_zoom_vals(
            vals,
            self.options,
            BigBedWrite::process_chrom_zoom,
            &runtime,
            &chrom_ids,
            (summary.bases_covered as f64 / summary.total_items as f64) as u32,
            zoom_counts,
            file,
            data_size,
        );
        let (mut file, zoom_entries, zoom_uncompress_buf_size) = output?;
        uncompress_buf_size = uncompress_buf_size.max(zoom_uncompress_buf_size);
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

    async fn process_val<I: ChromValues<Value = BedEntry>>(
        current_val: BedEntry,
        chrom_length: u32,
        chrom: &String,
        chrom_values: &mut I,
        summary: &mut Option<Summary>,
        state_val: &mut EntriesSection,
        options: BBIWriteOptions,
        runtime: &Handle,
        ftx: &mut ChromProcessingInputSectionChannel,
        chrom_id: u32,
    ) -> Result<(), ProcessChromError<I::Error>> {
        // Check a few preconditions:
        // - The current end is greater than or equal to the start
        // - The current end is at most the chromosome length
        // - If there is a next value, then it does not overlap value
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
        match chrom_values.peek() {
            None | Some(Err(_)) => (),
            Some(Ok(next_val)) => {
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
        }

        // Now, actually process the value.

        // First, update the summary.
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
            summary,
            current_val.start,
            current_val.end,
            chrom_values.peek().and_then(|v| v.ok()).map(|v| v.start),
        );

        // Then, add the current item to the actual values, and encode if full, or last item
        state_val.items.push(current_val);
        if chrom_values.peek().is_none() || state_val.items.len() >= options.items_per_slot as usize
        {
            let items = std::mem::replace(
                &mut state_val.items,
                Vec::with_capacity(options.items_per_slot as usize),
            );
            let handle = runtime
                .spawn(encode_section(options.compress, items, chrom_id))
                .map(|f| f.unwrap());
            ftx.send(handle.boxed()).await.expect("Couldn't send");
        }

        Ok(())
    }

    async fn process_val_zoom<I: ChromValues<Value = BedEntry>>(
        zoom_items: &mut Vec<ZoomItem>,
        options: BBIWriteOptions,
        item_start: u32,
        item_end: u32,
        chrom_values: &mut I,
        runtime: &Handle,
        chrom_id: u32,
    ) -> Result<(), ProcessChromError<I::Error>> {
        // Then, add the item to the zoom item queues. This is a bit complicated.
        for zoom_item in zoom_items.iter_mut() {
            debug_assert_ne!(zoom_item.records.len(), options.items_per_slot as usize);

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

            let next_start = chrom_values
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
                        if chrom_values.peek().is_none() {
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

        Ok(())
    }

    async fn process_chrom<I>(
        zooms_channels: Vec<(u32, ChromProcessingInputSectionChannel)>,
        mut ftx: ChromProcessingInputSectionChannel,
        chrom_id: u32,
        options: BBIWriteOptions,
        runtime: Handle,
        mut chrom_values: I,
        chrom: String,
        chrom_length: u32,
    ) -> Result<Summary, ProcessChromError<I::Error>>
    where
        I: ChromValues<Value = BedEntry>,
    {
        let mut summary: Option<Summary> = None;

        let mut state_val = EntriesSection {
            items: Vec::with_capacity(options.items_per_slot as usize),
            overlap: IndexList::new(),
        };
        let mut zoom_items = zooms_channels
            .into_iter()
            .map(|(size, channel)| ZoomItem {
                size,
                live_info: None,
                overlap: IndexList::new(),
                records: Vec::with_capacity(options.items_per_slot as usize),
                channel,
            })
            .collect();
        let mut total_items = 0;
        while let Some(current_val) = chrom_values.next() {
            // If there is a source error, propogate that up
            let current_val = current_val.map_err(ProcessChromError::SourceError)?;
            total_items += 1;

            let item_start = current_val.start;
            let item_end = current_val.end;

            BigBedWrite::process_val(
                current_val,
                chrom_length,
                &chrom,
                &mut chrom_values,
                &mut summary,
                &mut state_val,
                options,
                &runtime,
                &mut ftx,
                chrom_id,
            )
            .await?;

            BigBedWrite::process_val_zoom(
                &mut zoom_items,
                options,
                item_start,
                item_end,
                &mut chrom_values,
                &runtime,
                chrom_id,
            )
            .await?;
        }

        debug_assert!(state_val.items.is_empty());
        for zoom_item in zoom_items.iter_mut() {
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

    pub(crate) async fn process_chrom_no_zooms<I: ChromValues<Value = BedEntry>>(
        mut ftx: ChromProcessingInputSectionChannel,
        chrom_id: u32,
        options: BBIWriteOptions,
        runtime: Handle,
        mut chrom_values: I,
        chrom: String,
        chrom_length: u32,
    ) -> Result<(Summary, Vec<(u64, u64)>), ProcessChromError<I::Error>> {
        #[derive(Debug, Copy, Clone)]
        struct ZoomCounts {
            resolution: u64,
            current_end: u64,
            counts: u64,
        }

        let mut summary: Option<Summary> = None;

        let mut state_val = EntriesSection {
            items: Vec::with_capacity(options.items_per_slot as usize),
            overlap: IndexList::new(),
        };
        let mut zoom_counts: Vec<ZoomCounts> = std::iter::successors(Some(10), |z| Some(z * 4))
            .take_while(|z| *z <= u64::MAX / 4 && *z <= chrom_length as u64 * 4)
            .map(|z| ZoomCounts {
                resolution: z,
                current_end: 0,
                counts: 0,
            })
            .collect();

        let mut total_items = 0;
        while let Some(current_val) = chrom_values.next() {
            // If there is a source error, propogate that up
            let current_val = current_val.map_err(ProcessChromError::SourceError)?;
            total_items += 1;

            let item_start = current_val.start;
            let item_end = current_val.end;

            BigBedWrite::process_val(
                current_val,
                chrom_length,
                &chrom,
                &mut chrom_values,
                &mut summary,
                &mut state_val,
                options,
                &runtime,
                &mut ftx,
                chrom_id,
            )
            .await?;

            for zoom in &mut zoom_counts {
                if item_start as u64 >= zoom.current_end {
                    zoom.counts += 1;
                    zoom.current_end = item_start as u64 + zoom.resolution;
                }
                while item_end as u64 > zoom.current_end {
                    zoom.counts += 1;
                    zoom.current_end += zoom.resolution;
                }
            }
        }

        debug_assert!(state_val.items.is_empty());

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

        let zoom_counts = zoom_counts
            .into_iter()
            .map(|z| (z.resolution, z.counts))
            .collect();

        Ok((summary_complete, zoom_counts))
    }

    pub(crate) async fn process_chrom_zoom<I: ChromValues<Value = BedEntry>>(
        zooms_channels: Vec<(u32, ChromProcessingInputSectionChannel)>,
        chrom_id: u32,
        options: BBIWriteOptions,
        runtime: Handle,
        mut chrom_values: I,
    ) -> Result<(), ProcessChromError<I::Error>> {
        let mut zoom_items: Vec<ZoomItem> = zooms_channels
            .into_iter()
            .map(|(size, channel)| ZoomItem {
                size,
                live_info: None,
                overlap: IndexList::new(),
                records: Vec::with_capacity(options.items_per_slot as usize),
                channel,
            })
            .collect();

        while let Some(current_val) = chrom_values.next() {
            // If there is a source error, propogate that up
            let current_val = current_val.map_err(ProcessChromError::SourceError)?;

            let item_start = current_val.start;
            let item_end = current_val.end;

            BigBedWrite::process_val_zoom(
                &mut zoom_items,
                options,
                item_start,
                item_end,
                &mut chrom_values,
                &runtime,
                chrom_id,
            )
            .await?;
        }

        for zoom_item in zoom_items.iter_mut() {
            debug_assert!(zoom_item.live_info.is_none());
            debug_assert!(zoom_item.records.is_empty());
        }

        Ok(())
    }
}

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
