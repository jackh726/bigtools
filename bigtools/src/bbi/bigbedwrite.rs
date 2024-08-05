use std::collections::HashMap;
use std::ffi::CString;
use std::fs::File;
use std::io::{self, BufWriter, Seek, Write};
use std::path::Path;

use futures::sink::SinkExt;

use byteorder::{NativeEndian, WriteBytesExt};
use tokio::runtime::{Handle, Runtime};

use crate::bbiwrite::process_internal::BBIDataProcessorCreate;
use crate::utils::tell::Tell;
use crate::{
    write_info, BBIDataProcessor, BBIDataProcessoredData, BBIDataProcessoringInputSectionChannel,
    BBIDataSource, InternalProcessData, InternalTempZoomInfo, NoZoomsInternalProcessData,
    NoZoomsInternalProcessedData, ProcessDataError, ZoomsInternalProcessData,
    ZoomsInternalProcessedData,
};
use index_list::IndexList;

use crate::bbi::{BedEntry, Summary, Value, ZoomRecord, BIGBED_MAGIC};
use crate::bbiwrite::{
    self, encode_zoom_section, write_blank_headers, write_zooms, BBIProcessError, BBIWriteOptions,
    SectionData,
};
use crate::bed::autosql::parse::parse_autosql;

/// The struct used to write a bigBed file
pub struct BigBedWrite<W: Write + Seek + Send + 'static> {
    out: W,
    chrom_sizes: HashMap<String, u32>,
    pub options: BBIWriteOptions,
    pub autosql: Option<String>,
}

impl BigBedWrite<File> {
    pub fn create_file(
        path: impl AsRef<Path>,
        chrom_sizes: HashMap<String, u32>,
    ) -> io::Result<Self> {
        let out = File::create(path)?;
        Ok(BigBedWrite::new(out, chrom_sizes))
    }
}

impl<W: Write + Seek + Send + 'static> BigBedWrite<W> {
    pub fn new(out: W, chrom_sizes: HashMap<String, u32>) -> Self {
        BigBedWrite {
            out,
            chrom_sizes,
            options: BBIWriteOptions::default(),
            autosql: None,
        }
    }

    fn write_pre(
        file: &mut BufWriter<W>,
        autosql: Option<String>,
    ) -> Result<(u64, u64, u64, u64, u16), ProcessDataError> {
        write_blank_headers(file)?;

        let autosql = autosql.unwrap_or_else(|| crate::bed::autosql::BED3.to_string());

        let field_count = 'field_count: {
            let Ok(mut declarations) = parse_autosql(&autosql) else {
                break 'field_count None;
            };
            let Some(decl) = declarations.pop() else {
                break 'field_count None;
            };
            Some(decl.fields.len())
        };
        let field_count = field_count.unwrap_or(3) as u16;

        let autosql = CString::new(autosql.into_bytes()).map_err(|_| {
            ProcessDataError::InvalidInput("Invalid autosql: null byte in string".to_owned())
        })?;

        let autosql_offset = file.tell()?;
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
            field_count,
        ))
    }

    /// Write the values from `V` as a bigWig. Will utilize the provided runtime for encoding values and for reading through the values (potentially parallelized by chromosome).
    pub fn write<V: BBIDataSource<Value = BedEntry>>(
        self,
        vals: V,
        runtime: Runtime,
    ) -> Result<(), BBIProcessError<V::Error>> {
        let mut file = BufWriter::new(self.out);

        let (autosql_offset, total_summary_offset, full_data_offset, pre_data, field_count) =
            BigBedWrite::write_pre(&mut file, self.autosql)?;

        let output = bbiwrite::write_vals::<_, _, BigBedFullProcess>(
            vals,
            file,
            self.options,
            runtime,
            &self.chrom_sizes,
        );
        let (chrom_ids, summary, mut file, raw_sections_iter, zoom_infos, uncompress_buf_size) =
            output?;

        let chrom_ids = chrom_ids.get_map();
        let (data_size, chrom_index_start, index_start, _total_sections) = bbiwrite::write_mid(
            &mut file,
            pre_data,
            raw_sections_iter,
            self.chrom_sizes,
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
            field_count,
            // No separate option to specify field count, so use field count
            // defined in autosql
            field_count,
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
    pub fn write_multipass<V: BBIDataSource<Value = BedEntry>>(
        self,
        make_vals: impl Fn() -> Result<V, BBIProcessError<V::Error>>,
        runtime: Runtime,
    ) -> Result<(), BBIProcessError<V::Error>> {
        let mut file = BufWriter::new(self.out);

        let (autosql_offset, total_summary_offset, full_data_offset, pre_data, field_count) =
            BigBedWrite::write_pre(&mut file, self.autosql)?;

        let vals = make_vals()?;

        let output = bbiwrite::write_vals_no_zoom::<_, _, BigBedNoZoomsProcess>(
            vals,
            file,
            self.options,
            &runtime,
            &self.chrom_sizes,
        );
        let (chrom_ids, summary, zoom_counts, mut file, raw_sections_iter, mut uncompress_buf_size) =
            output?;

        let chrom_ids = chrom_ids.get_map();
        let (data_size, chrom_index_start, index_start, _total_sections) = bbiwrite::write_mid(
            &mut file,
            pre_data,
            raw_sections_iter,
            self.chrom_sizes,
            &chrom_ids,
            self.options,
        )?;

        let vals = make_vals()?;

        let output = bbiwrite::write_zoom_vals::<_, _, BigBedZoomsProcess<W>>(
            vals,
            self.options,
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
            field_count,
            // No separate option to specify field count, so use field count
            // defined in autosql
            field_count,
            autosql_offset,
            total_summary_offset,
            uncompress_buf_size,
            zoom_entries,
            summary,
            summary.total_items,
        )?;

        Ok(())
    }
}

async fn process_val(
    current_val: BedEntry,
    next_val: Option<&BedEntry>,
    chrom_length: u32,
    chrom: &String,
    summary: &mut Option<Summary>,
    items: &mut Vec<BedEntry>,
    overlap: &mut IndexList<Value>,
    options: BBIWriteOptions,
    runtime: &Handle,
    ftx: &mut BBIDataProcessoringInputSectionChannel,
    chrom_id: u32,
) -> Result<(), ProcessDataError> {
    // Check a few preconditions:
    // - The current end is greater than or equal to the start
    // - The current end is at most the chromosome length
    // - If there is a next value, then it does not overlap value
    // TODO: test these correctly fails
    if current_val.start > current_val.end {
        return Err(ProcessDataError::InvalidInput(format!(
            "Invalid bed: {} > {}",
            current_val.start, current_val.end
        )));
    }
    if current_val.start >= chrom_length {
        return Err(ProcessDataError::InvalidInput(format!(
            "Invalid bed: `{}` is greater than the chromosome ({}) length ({})",
            current_val.start, chrom, chrom_length
        )));
    }
    match next_val {
        None => (),
        Some(next_val) => {
            if current_val.start > next_val.start {
                return Err(ProcessDataError::InvalidInput(format!(
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
                .get_first()
                .map(|f| f.start == item_start)
                .unwrap_or(true));

            // For each item in `overlap` that overlaps the current
            // item, add `1` to the value.
            let mut index = overlap.first_index();
            while index.is_some() {
                match overlap.get_mut(index) {
                    None => break,
                    Some(o) => {
                        o.value += 1.0;
                        if item_end < o.end {
                            let value = o.value - 1.0;
                            let end = o.end;
                            o.end = item_end;
                            overlap.insert_after(
                                index,
                                Value {
                                    start: item_end,
                                    end,
                                    value,
                                },
                            );
                            break;
                        }
                        index = overlap.next_index(index);
                    }
                }
            }

            debug_assert!(overlap
                .get_last()
                .map(|o| o.end >= item_start)
                .unwrap_or(true));

            if overlap.get_last().map(|o| o.end).unwrap_or(item_start) == item_start {
                overlap.insert_last(Value {
                    start: item_start,
                    end: item_end,
                    value: 1.0,
                });
            }

            let next_start = next_start_opt.unwrap_or(u32::max_value());

            while overlap
                .get_first()
                .map(|f| f.start < next_start)
                .unwrap_or(false)
            {
                let mut removed = overlap.remove_first().unwrap();
                let (len, val) = if removed.end <= next_start {
                    (removed.end - removed.start, f64::from(removed.value))
                } else {
                    let len = next_start - removed.start;
                    let val = f64::from(removed.value);
                    removed.start = next_start;
                    overlap.insert_first(removed);
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
        overlap,
        summary,
        current_val.start,
        current_val.end,
        next_val.map(|v| v.start),
    );

    // Then, add the current item to the actual values, and encode if full, or last item
    items.push(current_val);
    if next_val.is_none() || items.len() >= options.items_per_slot as usize {
        let items = std::mem::replace(items, Vec::with_capacity(options.items_per_slot as usize));
        let handle = runtime.spawn(encode_section(options.compress, items, chrom_id));
        ftx.send(handle).await.expect("Couldn't send");
    }

    Ok(())
}

async fn process_val_zoom(
    zoom_items: &mut Vec<ZoomItem>,
    options: BBIWriteOptions,
    item_start: u32,
    item_end: u32,
    next_val: Option<&BedEntry>,
    runtime: &Handle,
    chrom_id: u32,
) -> Result<(), ProcessDataError> {
    // Then, add the item to the zoom item queues. This is a bit complicated.
    for zoom_item in zoom_items.iter_mut() {
        debug_assert_ne!(zoom_item.records.len(), options.items_per_slot as usize);

        let overlap = &mut zoom_item.overlap;

        // For each item in `overlap` that overlaps the current
        // item, add `1` to the value.
        let mut index = overlap.first_index();
        while index.is_some() {
            match overlap.get_mut(index) {
                None => break,
                Some(o) => {
                    o.value += 1.0;
                    if item_end < o.end {
                        let value = o.value - 1.0;
                        let end = o.end;
                        o.end = item_end;
                        overlap.insert_after(
                            index,
                            Value {
                                start: item_end,
                                end,
                                value,
                            },
                        );
                        break;
                    }
                    index = overlap.next_index(index);
                }
            }
        }

        debug_assert!(overlap
            .get_last()
            .map(|o| o.end >= item_start)
            .unwrap_or(true));

        if overlap.get_last().map(|o| o.end).unwrap_or(item_start) == item_start {
            overlap.insert_last(Value {
                start: item_start,
                end: item_end,
                value: 1.0,
            });
        }

        let next_start = next_val.map(|v| v.start).unwrap_or(u32::max_value());

        while overlap
            .get_first()
            .map(|f| f.start < next_start)
            .unwrap_or(false)
        {
            let mut removed = overlap.remove_first().unwrap();
            let val = f64::from(removed.value);
            let (removed_start, removed_end) = if removed.end <= next_start {
                (removed.start, removed.end)
            } else {
                let start = removed.start;
                removed.start = next_start;
                overlap.insert_first(removed);
                (start, next_start)
            };

            let mut add_start = removed_start;
            loop {
                if add_start >= removed_end {
                    if next_val.is_none() {
                        if let Some((mut zoom2, total_items)) = zoom_item.live_info.take() {
                            zoom2.summary.total_items = total_items;
                            zoom_item.records.push(zoom2);
                        }
                        if !zoom_item.records.is_empty() {
                            let items = std::mem::take(&mut zoom_item.records);
                            let handle =
                                runtime.spawn(encode_zoom_section(options.compress, items));
                            zoom_item.channel.send(handle).await.expect("Couln't send");
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
                    let handle = runtime.spawn(encode_zoom_section(options.compress, items));
                    zoom_item.channel.send(handle).await.expect("Couln't send");
                }
            }
        }

        debug_assert_ne!(zoom_item.records.len(), options.items_per_slot as usize);
    }

    Ok(())
}

// While we do technically lose precision here by using the f32 in Value, we can reuse the same merge_into method
struct ZoomItem {
    size: u32,
    live_info: Option<(ZoomRecord, u64)>,
    overlap: IndexList<Value>,
    records: Vec<ZoomRecord>,
    channel: BBIDataProcessoringInputSectionChannel,
}
struct EntriesSection {
    items: Vec<BedEntry>,
    overlap: IndexList<Value>,
    zoom_items: Vec<ZoomItem>,
}

pub(crate) struct BigBedFullProcess {
    summary: Option<Summary>,
    state_val: EntriesSection,
    total_items: u64,

    ftx: BBIDataProcessoringInputSectionChannel,
    chrom_id: u32,
    options: BBIWriteOptions,
    runtime: Handle,
    chrom: String,
    length: u32,
}

impl BBIDataProcessorCreate for BigBedFullProcess {
    type I = InternalProcessData;
    type Out = BBIDataProcessoredData;
    fn destroy(self) -> BBIDataProcessoredData {
        let Self {
            summary,
            total_items,
            state_val,
            ..
        } = self;

        debug_assert!(state_val.items.is_empty());
        for zoom_item in state_val.zoom_items.iter() {
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
        BBIDataProcessoredData(summary_complete)
    }
    fn create(internal_data: InternalProcessData) -> Self {
        let InternalProcessData(zooms_channels, ftx, chrom_id, options, runtime, chrom, length) =
            internal_data;

        let summary: Option<Summary> = None;

        let zoom_items = zooms_channels
            .into_iter()
            .map(|(size, channel)| ZoomItem {
                size,
                live_info: None,
                overlap: IndexList::new(),
                records: Vec::with_capacity(options.items_per_slot as usize),
                channel,
            })
            .collect();
        let state_val = EntriesSection {
            zoom_items,
            items: Vec::with_capacity(options.items_per_slot as usize),
            overlap: IndexList::new(),
        };
        let total_items = 0;
        BigBedFullProcess {
            summary,
            state_val,
            total_items,
            ftx,
            chrom_id,
            options,
            runtime,
            chrom,
            length,
        }
    }
}
impl BBIDataProcessor for BigBedFullProcess {
    type Value = BedEntry;
    async fn do_process(
        &mut self,
        current_val: Self::Value,
        next_val: Option<&Self::Value>,
    ) -> Result<(), ProcessDataError> {
        let Self {
            summary,
            total_items,
            state_val,
            ftx,
            chrom_id,
            options,
            runtime,
            chrom,
            length,
        } = self;
        let chrom_id = *chrom_id;
        let length = *length;

        *total_items += 1;

        let item_start = current_val.start;
        let item_end = current_val.end;

        process_val(
            current_val,
            next_val,
            length,
            &chrom,
            summary,
            &mut state_val.items,
            &mut state_val.overlap,
            *options,
            &runtime,
            ftx,
            chrom_id,
        )
        .await?;

        process_val_zoom(
            &mut state_val.zoom_items,
            *options,
            item_start,
            item_end,
            next_val,
            &runtime,
            chrom_id,
        )
        .await?;

        Ok(())
    }
}

#[derive(Debug, Copy, Clone)]
struct ZoomCounts {
    resolution: u64,
    current_end: u64,
    counts: u64,
}
struct BigBedNoZoomsProcess {
    ftx: BBIDataProcessoringInputSectionChannel,
    chrom_id: u32,
    options: BBIWriteOptions,
    runtime: Handle,
    chrom: String,
    length: u32,

    summary: Option<Summary>,
    items: Vec<BedEntry>,
    overlap: IndexList<Value>,
    zoom_counts: Vec<ZoomCounts>,
    total_items: u64,
}

impl BBIDataProcessorCreate for BigBedNoZoomsProcess {
    type I = NoZoomsInternalProcessData;
    type Out = NoZoomsInternalProcessedData;
    fn create(internal_data: Self::I) -> Self {
        let NoZoomsInternalProcessData(ftx, chrom_id, options, runtime, chrom, length) =
            internal_data;

        let summary = None;

        let items: Vec<BedEntry> = Vec::with_capacity(options.items_per_slot as usize);
        let zoom_counts: Vec<ZoomCounts> = std::iter::successors(Some(10), |z| Some(z * 4))
            .take_while(|z| *z <= u64::MAX / 4 && *z <= length as u64 * 4)
            .map(|z| ZoomCounts {
                resolution: z,
                current_end: 0,
                counts: 0,
            })
            .collect();

        BigBedNoZoomsProcess {
            ftx,
            chrom_id,
            options,
            runtime,
            chrom,
            length,
            summary,
            items,
            overlap: IndexList::new(),
            zoom_counts,
            total_items: 0,
        }
    }
    fn destroy(self) -> Self::Out {
        let BigBedNoZoomsProcess {
            items,
            summary,
            zoom_counts,
            total_items,
            ..
        } = self;

        debug_assert!(items.is_empty());

        let mut summary = summary.unwrap_or(Summary {
            total_items: 0,
            bases_covered: 0,
            min_val: 0.0,
            max_val: 0.0,
            sum: 0.0,
            sum_squares: 0.0,
        });
        summary.total_items = total_items;

        let zoom_counts = zoom_counts
            .into_iter()
            .map(|z| (z.resolution, z.counts))
            .collect();

        NoZoomsInternalProcessedData(summary, zoom_counts)
    }
}

impl BBIDataProcessor for BigBedNoZoomsProcess {
    type Value = BedEntry;
    async fn do_process(
        &mut self,
        current_val: Self::Value,
        next_val: Option<&Self::Value>,
    ) -> Result<(), ProcessDataError> {
        let BigBedNoZoomsProcess {
            ftx,
            chrom_id,
            options,
            runtime,
            chrom,
            length,
            summary,
            items,
            overlap,
            zoom_counts,
            total_items,
        } = self;

        *total_items += 1;

        let item_start = current_val.start;
        let item_end = current_val.end;

        process_val(
            current_val,
            next_val,
            *length,
            &chrom,
            summary,
            items,
            overlap,
            *options,
            &runtime,
            ftx,
            *chrom_id,
        )
        .await?;

        for zoom in zoom_counts {
            if item_start as u64 >= zoom.current_end {
                zoom.counts += 1;
                zoom.current_end = item_start as u64 + zoom.resolution;
            }
            while item_end as u64 > zoom.current_end {
                zoom.counts += 1;
                zoom.current_end += zoom.resolution;
            }
        }

        Ok(())
    }
}

struct BigBedZoomsProcess<W: Write + Seek + Send + 'static> {
    temp_zoom_items: Vec<InternalTempZoomInfo<W>>,
    chrom_id: u32,
    options: BBIWriteOptions,
    runtime: Handle,

    zoom_items: Vec<ZoomItem>,
}

impl<W: Write + Seek + Send + 'static> BBIDataProcessorCreate for BigBedZoomsProcess<W> {
    type I = ZoomsInternalProcessData<W>;
    type Out = ZoomsInternalProcessedData<W>;
    fn create(internal_data: Self::I) -> Self {
        let ZoomsInternalProcessData(temp_zoom_items, zooms_channels, chrom_id, options, runtime) =
            internal_data;

        let zoom_items: Vec<ZoomItem> = zooms_channels
            .into_iter()
            .map(|(size, channel)| ZoomItem {
                size,
                live_info: None,
                overlap: IndexList::new(),
                records: Vec::with_capacity(options.items_per_slot as usize),
                channel,
            })
            .collect();

        BigBedZoomsProcess {
            temp_zoom_items,
            chrom_id,
            options,
            runtime,
            zoom_items,
        }
    }
    fn destroy(self) -> Self::Out {
        let BigBedZoomsProcess { zoom_items, .. } = self;

        for zoom_item in zoom_items.iter() {
            debug_assert!(zoom_item.live_info.is_none());
            debug_assert!(zoom_item.records.is_empty());
        }

        ZoomsInternalProcessedData(self.temp_zoom_items)
    }
}
impl<W: Write + Seek + Send + 'static> BBIDataProcessor for BigBedZoomsProcess<W> {
    type Value = BedEntry;
    async fn do_process(
        &mut self,
        current_val: Self::Value,
        next_val: Option<&Self::Value>,
    ) -> Result<(), ProcessDataError> {
        let BigBedZoomsProcess {
            chrom_id,
            options,
            runtime,
            zoom_items,
            ..
        } = self;

        process_val_zoom(
            zoom_items,
            *options,
            current_val.start,
            current_val.end,
            next_val,
            &runtime,
            *chrom_id,
        )
        .await?;

        Ok(())
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
