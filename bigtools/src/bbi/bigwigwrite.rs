/*!
Provides the interface for writing bigWig files.

## Example
```rust,no_run
# use std::collections::HashMap;
# use std::error::Error;
# use std::path::PathBuf;
# use std::fs::File;
# use bigtools::BigWigWrite;
# use bigtools::bedchromdata::BedParserStreamingIterator;
# use bigtools::bed::bedparser::BedParser;
# fn main() -> Result<(), Box<dyn Error>> {
# let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
# dir.push("resources/test");
# let mut bedgraph_in = dir.clone();
# bedgraph_in.push("single_chrom.bedGraph");
// First, set up our input data. Here, we're using the `BedParserStreamingIterator` with a `BedParser`.
let bedgraph_file: File = File::open(bedgraph_in)?;
let vals_iter = BedParser::from_bedgraph_file(bedgraph_file);
let vals = BedParserStreamingIterator::new(vals_iter, false);

// Then, we need to know what the chromosome sizes are. This can be read in from a file, but here we
// just construct a map for ease.
let mut chrom_map = HashMap::new();
chrom_map.insert("chr17".to_string(), 83257441);

// We also need a `Runtime` to spawn processing on.
let runtime = tokio::runtime::Builder::new_multi_thread()
    .worker_threads(6)
    .build()
    .expect("Unable to create runtime.");

// Finally, we can create a `BigWigWrite` with a file to write to. We'll use a temporary file.
let tempfile = tempfile::NamedTempFile::new()?;
let out = BigWigWrite::create_file(tempfile.path().to_string_lossy().to_string());
// Then write.
out.write(chrom_map, vals, runtime)?;
# Ok(())
# }
```
*/
use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fs::File;
use std::future::Future;
use std::io::{self, BufWriter, Write};
use std::vec;

use futures::channel::mpsc as futures_mpsc;
use futures::future::FutureExt;
use futures::sink::SinkExt;

use byteorder::{NativeEndian, WriteBytesExt};
use tokio::runtime::{Handle, Runtime};

use crate::bbiwrite::process_internal::ChromProcessCreate;
use crate::utils::chromvalues::ChromValues;
use crate::utils::idmap::IdMap;
use crate::utils::tell::Tell;
use crate::utils::tempfilebuffer::{TempFileBuffer, TempFileBufferWriter};
use crate::{
    future_channel, write_chroms_with_zooms, write_info, ChromData, ChromData2, ChromProcess,
    ChromProcessedData, ChromProcessingInputSectionChannel, InternalProcessData, Section,
    TempZoomInfo, ZoomInfo, ZoomValue,
};

use crate::bbi::{Summary, Value, ZoomRecord, BIGWIG_MAGIC};
use crate::bbiwrite::{
    self, encode_zoom_section, get_rtreeindex, write_blank_headers, write_chrom_tree,
    write_rtreeindex, write_zooms, BBIWriteOptions, ProcessChromError, SectionData,
};

struct ZoomItem {
    // How many bases this zoom item covers
    size: u32,
    // The current zoom entry
    live_info: Option<ZoomRecord>,
    // All zoom entries in the current section
    records: Vec<ZoomRecord>,
    channel: ChromProcessingInputSectionChannel,
}

/// The struct used to write a bigWig file
pub struct BigWigWrite {
    pub path: String,
    pub options: BBIWriteOptions,
}

impl BigWigWrite {
    pub fn create_file(path: String) -> Self {
        BigWigWrite {
            path,
            options: BBIWriteOptions::default(),
        }
    }

    fn write_pre<E: Error>(
        file: &mut BufWriter<File>,
    ) -> Result<(u64, u64, u64), ProcessChromError<E>> {
        write_blank_headers(file)?;

        let total_summary_offset = file.tell()?;
        file.write_all(&[0; 40])?;

        let full_data_offset = file.tell()?;

        // Total items
        // Unless we know the vals ahead of time, we can't estimate total sections ahead of time.
        // Even then simply doing "(vals.len() as u32 + ITEMS_PER_SLOT - 1) / ITEMS_PER_SLOT"
        // underestimates because sections are split by chrom too, not just size.
        // Skip for now, and come back when we write real header + summary.
        file.write_u64::<NativeEndian>(0)?;

        let pre_data = file.tell()?;

        Ok((total_summary_offset, full_data_offset, pre_data))
    }

    fn write_mid<E: Error>(
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

        // Since the chrom tree is read before the index, we put this before the full data index
        // Therefore, there is a higher likelihood that the udc file will only need one read for chrom tree + full data index
        // Putting the chrom tree before the data also has a higher likelihood of being included with the beginning headers,
        // but requires us to know all the data ahead of time (when writing)
        let chrom_index_start = file.tell()?;
        write_chrom_tree(file, chrom_sizes, chrom_ids)?;

        let index_start = file.tell()?;
        let (nodes, levels, total_sections) = get_rtreeindex(sections_iter, options);
        write_rtreeindex(file, nodes, levels, total_sections, options)?;

        Ok((data_size, chrom_index_start, index_start, total_sections))
    }

    /// Write the values from `V` as a bigWig. Will utilize the provided runtime for encoding values and for reading through the values (potentially parallelized by chromosome).
    pub fn write<
        Values: ChromValues<Value = Value> + Send + 'static,
        V: ChromData2<Values = Values>,
    >(
        self,
        chrom_sizes: HashMap<String, u32>,
        mut vals: V,
        runtime: Runtime,
    ) -> Result<(), ProcessChromError<Values::Error>> {
        let options = self.options;
        let fp = File::create(self.path.clone())?;
        let mut file = BufWriter::new(fp);

        let (total_summary_offset, full_data_offset, pre_data) = BigWigWrite::write_pre(&mut file)?;

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

        let mut summary: Option<Summary> = None;
        let (send, recv) = futures_mpsc::unbounded();
        let write_fut = write_chroms_with_zooms(file, zooms_map, recv);
        let (write_fut, write_fut_handle) = write_fut.remote_handle();
        runtime.spawn(write_fut);

        let handle = runtime.handle();

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
        let mut do_read = |chrom: String| -> Result<_, ProcessChromError<_>> {
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

            let internal_data = crate::InternalProcessData(
                zooms_channels,
                ftx,
                chrom_id,
                options,
                runtime.handle().clone(),
                chrom,
                length,
            );
            Ok(BigWigFullProcess::create(internal_data))
        };

        let mut advance = |p: BigWigFullProcess| {
            let data = p.destroy();
            let ChromProcessedData(chrom_summary) = data;
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
        };

        vals.process_to_bbi(handle, &mut do_read, &mut advance)?;

        drop(send);

        self.write_internal_post(
            summary,
            runtime,
            write_fut_handle,
            chrom_ids,
            pre_data,
            chrom_sizes,
            full_data_offset,
            total_summary_offset,
        )
    }

    /// Write the values from `V` as a bigWig. Will utilize the provided runtime for encoding values, but will read through values on the current thread.
    pub fn write_singlethreaded<
        Values: ChromValues<Value = Value>,
        V: ChromData2<Values = Values>,
    >(
        self,
        chrom_sizes: HashMap<String, u32>,
        mut vals: V,
        runtime: Runtime,
    ) -> Result<(), ProcessChromError<Values::Error>> {
        let options = self.options;
        let fp = File::create(self.path.clone())?;
        let mut file = BufWriter::new(fp);

        let (total_summary_offset, full_data_offset, pre_data) = BigWigWrite::write_pre(&mut file)?;

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

        let mut summary: Option<Summary> = None;
        let (send, recv) = futures_mpsc::unbounded();
        let write_fut = write_chroms_with_zooms(file, zooms_map, recv);
        let (write_fut, write_fut_handle) = write_fut.remote_handle();
        runtime.spawn(write_fut);

        let handle = runtime.handle();

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
        let mut do_read = |chrom: String| -> Result<_, ProcessChromError<_>> {
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

            let internal_data = crate::InternalProcessData(
                zooms_channels,
                ftx,
                chrom_id,
                options,
                runtime.handle().clone(),
                chrom,
                length,
            );
            Ok(BigWigFullProcess::create(internal_data))
        };

        let mut advance = |p: BigWigFullProcess| {
            let data = p.destroy();
            let ChromProcessedData(chrom_summary) = data;
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
        };

        vals.process_to_bbi(handle, &mut do_read, &mut advance)?;

        drop(send);

        self.write_internal_post(
            summary,
            runtime,
            write_fut_handle,
            chrom_ids,
            pre_data,
            chrom_sizes,
            full_data_offset,
            total_summary_offset,
        )
    }

    fn write_internal_post<E: Error>(
        self,
        summary: Option<Summary>,
        runtime: Runtime,
        write_fut_handle: impl Future<
            Output = Result<
                (
                    BufWriter<File>,
                    usize,
                    Vec<crossbeam_channel::IntoIter<Section>>,
                    BTreeMap<u32, ZoomValue>,
                ),
                ProcessChromError<E>,
            >,
        >,
        chrom_ids: IdMap,
        pre_data: u64,
        chrom_sizes: HashMap<String, u32>,
        full_data_offset: u64,
        total_summary_offset: u64,
    ) -> Result<(), ProcessChromError<E>> {
        let summary = summary.unwrap_or(Summary {
            total_items: 0,
            bases_covered: 0,
            min_val: 0.0,
            max_val: 0.0,
            sum: 0.0,
            sum_squares: 0.0,
        });

        let (mut file, max_uncompressed_buf_size, section_iter, zooms_map) =
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
        let raw_sections_iter = section_iter.into_iter().flatten();

        let chrom_ids = chrom_ids.get_map();
        let (data_size, chrom_index_start, index_start, total_sections) = BigWigWrite::write_mid(
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
            BIGWIG_MAGIC,
            num_zooms,
            chrom_index_start,
            full_data_offset,
            index_start,
            0,
            0,
            0,
            total_summary_offset,
            max_uncompressed_buf_size,
            zoom_entries,
            summary,
            total_sections,
        )?;

        Ok(())
    }

    /// Write the values from `V` as a bigWig. Will utilize the provided runtime for encoding values and for reading through the values (potentially parallelized by chromosome).
    /// This will take two passes on the provided values: first to write the values themselves, then the zooms. This is beneficial over `write` on smaller files, where the encoding of
    /// high resolution zooms takes up a substantial portion of total processing time.
    pub fn write_multipass<
        Values: ChromValues<Value = Value> + Send + 'static,
        V: ChromData<Values = Values>,
    >(
        self,
        make_vals: impl Fn() -> Result<V, ProcessChromError<Values::Error>>,
        chrom_sizes: HashMap<String, u32>,
        runtime: Runtime,
    ) -> Result<(), ProcessChromError<Values::Error>> {
        let fp = File::create(self.path.clone())?;
        let mut file = BufWriter::new(fp);

        let (total_summary_offset, full_data_offset, pre_data) = BigWigWrite::write_pre(&mut file)?;

        let vals = make_vals()?;

        let runtime_handle = runtime.handle();

        let process_chrom = |ftx: ChromProcessingInputSectionChannel,
                             chrom_id: u32,
                             options: BBIWriteOptions,
                             runtime: Handle,
                             chrom_values: Values,
                             chrom: String,
                             chrom_length: u32| {
            let fut = BigWigWrite::process_chrom_no_zooms(
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
        let (data_size, chrom_index_start, index_start, total_sections) = BigWigWrite::write_mid(
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
            BigWigWrite::process_chrom_zoom,
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
            BIGWIG_MAGIC,
            num_zooms,
            chrom_index_start,
            full_data_offset,
            index_start,
            0,
            0,
            0,
            total_summary_offset,
            uncompress_buf_size,
            zoom_entries,
            summary,
            total_sections,
        )?;

        Ok(())
    }

    async fn process_val<E: Error>(
        current_val: Value,
        next_val: Option<&Value>,
        chrom_length: u32,
        chrom: &String,
        summary: &mut Summary,
        items: &mut Vec<Value>,
        options: BBIWriteOptions,
        runtime: &Handle,
        ftx: &mut ChromProcessingInputSectionChannel,
        chrom_id: u32,
    ) -> Result<(), ProcessChromError<E>> {
        // Check a few preconditions:
        // - The current end is greater than or equal to the start
        // - The current end is at most the chromosome length
        // - If there is a next value, then it does not overlap value
        // TODO: test these correctly fails
        if current_val.start > current_val.end {
            return Err(ProcessChromError::InvalidInput(format!(
                "Invalid bed graph: {} > {}",
                current_val.start, current_val.end
            )));
        }
        if current_val.end > chrom_length {
            return Err(ProcessChromError::InvalidInput(format!(
                "Invalid bed graph: `{}` is greater than the chromosome ({}) length ({})",
                current_val.end, chrom, chrom_length
            )));
        }
        match next_val {
            None => {}
            Some(next_val) => {
                if current_val.end > next_val.start {
                    return Err(ProcessChromError::InvalidInput(format!(
                        "Invalid bed graph: overlapping values on chromosome {} at {}-{} and {}-{}",
                        chrom, current_val.start, current_val.end, next_val.start, next_val.end,
                    )));
                }
            }
        }

        // Now, actually process the value.

        // First, update the summary.
        let len = current_val.end - current_val.start;
        let val = f64::from(current_val.value);
        summary.total_items += 1;
        summary.bases_covered += u64::from(len);
        summary.min_val = summary.min_val.min(val);
        summary.max_val = summary.max_val.max(val);
        summary.sum += f64::from(len) * val;
        summary.sum_squares += f64::from(len) * val * val;

        // Then, add the current item to the actual values, and encode if full, or last item
        items.push(current_val);
        if next_val.is_none() || items.len() >= options.items_per_slot as usize {
            let items = std::mem::take(items);
            let handle = runtime
                .spawn(encode_section(options.compress, items, chrom_id))
                .map(|f| f.unwrap());
            ftx.send(handle.boxed()).await.expect("Couldn't send");
        }

        Ok(())
    }

    async fn process_val_zoom<E: Error>(
        zoom_items: &mut Vec<ZoomItem>,
        options: BBIWriteOptions,
        current_val: Value,
        next_val: Option<&Value>,
        runtime: &Handle,
        chrom_id: u32,
    ) -> Result<(), ProcessChromError<E>> {
        // Then, add the item to the zoom item queues. This is a bit complicated.
        for zoom_item in zoom_items.iter_mut() {
            debug_assert_ne!(zoom_item.records.len(), options.items_per_slot as usize);

            // Zooms are comprised of a tiled set of summaries. Each summary spans a fixed length.
            // Zoom summaries are compressed similarly to main data, with a given items per slot.
            // It may be the case that our value spans across multiple zoom summaries, so this inner loop handles that.

            // `add_start` indicates where we are *currently* adding bases from (either the start of this item or in the middle, but beginning of another zoom section)
            let mut add_start = current_val.start;
            loop {
                // Write section if full; or if no next section, some items, and no current zoom record
                if (add_start >= current_val.end
                    && zoom_item.live_info.is_none()
                    && next_val.is_none()
                    && !zoom_item.records.is_empty())
                    || zoom_item.records.len() == options.items_per_slot as usize
                {
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
                if add_start >= current_val.end {
                    if next_val.is_none() {
                        if let Some(zoom2) = zoom_item.live_info.take() {
                            zoom_item.records.push(zoom2);
                            continue;
                        }
                    }
                    break;
                }
                let val = f64::from(current_val.value);
                let zoom2 = zoom_item.live_info.get_or_insert(ZoomRecord {
                    chrom: chrom_id,
                    start: add_start,
                    end: add_start,
                    summary: Summary {
                        total_items: 0,
                        bases_covered: 0,
                        min_val: val,
                        max_val: val,
                        sum: 0.0,
                        sum_squares: 0.0,
                    },
                });
                // The end of zoom record
                let next_end = zoom2.start + zoom_item.size;
                // End of bases that we could add
                let add_end = std::cmp::min(next_end, current_val.end);
                // If the last zoom ends before this value starts, we don't add anything
                if add_end >= add_start {
                    let added_bases = add_end - add_start;
                    zoom2.end = add_end;
                    zoom2.summary.total_items += 1;
                    zoom2.summary.bases_covered += u64::from(added_bases);
                    zoom2.summary.min_val = zoom2.summary.min_val.min(val);
                    zoom2.summary.max_val = zoom2.summary.max_val.max(val);
                    zoom2.summary.sum += f64::from(added_bases) * val;
                    zoom2.summary.sum_squares += f64::from(added_bases) * val * val;
                }
                // If we made it to the end of the zoom (whether it was because the zoom ended before this value started,
                // or we added to the end of the zoom), then write this zooms to the current section
                if add_end == next_end {
                    zoom_item.records.push(zoom_item.live_info.take().unwrap());
                }
                // Set where we would start for next time
                add_start = add_end;
            }
            debug_assert_ne!(zoom_item.records.len(), options.items_per_slot as usize);
        }

        Ok(())
    }

    pub(crate) async fn process_chrom_no_zooms<I: ChromValues<Value = Value>>(
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

        let mut summary = Summary {
            total_items: 0,
            bases_covered: 0,
            min_val: f64::MAX,
            max_val: f64::MIN,
            sum: 0.0,
            sum_squares: 0.0,
        };

        let mut items: Vec<Value> = Vec::with_capacity(options.items_per_slot as usize);
        let mut zoom_counts: Vec<ZoomCounts> = std::iter::successors(Some(10), |z| Some(z * 4))
            .take_while(|z| *z <= u64::MAX / 4 && *z <= chrom_length as u64 * 4)
            .map(|z| ZoomCounts {
                resolution: z,
                current_end: 0,
                counts: 0,
            })
            .collect();

        while let Some(current_val) = chrom_values.next() {
            // If there is a source error, propogate that up
            let current_val = current_val.map_err(ProcessChromError::SourceError)?;
            let next_val = match chrom_values.peek() {
                None | Some(Err(_)) => None,
                Some(Ok(v)) => Some(v),
            };

            BigWigWrite::process_val(
                current_val,
                next_val,
                chrom_length,
                &chrom,
                &mut summary,
                &mut items,
                options,
                &runtime,
                &mut ftx,
                chrom_id,
            )
            .await?;

            for zoom in &mut zoom_counts {
                if current_val.start as u64 >= zoom.current_end {
                    zoom.counts += 1;
                    zoom.current_end = current_val.start as u64 + zoom.resolution;
                }
                while current_val.end as u64 > zoom.current_end {
                    zoom.counts += 1;
                    zoom.current_end += zoom.resolution;
                }
            }
        }

        debug_assert!(items.is_empty());

        if summary.total_items == 0 {
            summary.min_val = 0.0;
            summary.max_val = 0.0;
        }

        let zoom_counts = zoom_counts
            .into_iter()
            .map(|z| (z.resolution, z.counts))
            .collect();

        Ok((summary, zoom_counts))
    }

    pub(crate) async fn process_chrom_zoom<I: ChromValues<Value = Value>>(
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
                records: Vec::with_capacity(options.items_per_slot as usize),
                channel,
            })
            .collect();

        while let Some(current_val) = chrom_values.next() {
            // If there is a source error, propogate that up
            let current_val = current_val.map_err(ProcessChromError::SourceError)?;
            let next_val = match chrom_values.peek() {
                None | Some(Err(_)) => None,
                Some(Ok(v)) => Some(v),
            };

            BigWigWrite::process_val_zoom(
                &mut zoom_items,
                options,
                current_val,
                next_val,
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

pub(crate) struct BigWigFullProcess {
    summary: Summary,
    items: Vec<Value>,
    zoom_items: Vec<ZoomItem>,

    ftx: ChromProcessingInputSectionChannel,
    chrom_id: u32,
    options: BBIWriteOptions,
    runtime: Handle,
    chrom: String,
    length: u32,
}

impl ChromProcessCreate for BigWigFullProcess {
    fn create(internal_data: InternalProcessData) -> Self {
        let InternalProcessData(zooms_channels, ftx, chrom_id, options, runtime, chrom, length) =
            internal_data;

        let summary = Summary {
            total_items: 0,
            bases_covered: 0,
            min_val: f64::MAX,
            max_val: f64::MIN,
            sum: 0.0,
            sum_squares: 0.0,
        };

        let items = Vec::with_capacity(options.items_per_slot as usize);
        let zoom_items: Vec<ZoomItem> = zooms_channels
            .into_iter()
            .map(|(size, channel)| ZoomItem {
                size,
                live_info: None,
                records: Vec::with_capacity(options.items_per_slot as usize),
                channel,
            })
            .collect();

        BigWigFullProcess {
            summary,
            items,
            zoom_items,
            ftx,
            chrom_id,
            options,
            runtime,
            chrom,
            length,
        }
    }
    fn destroy(self) -> ChromProcessedData {
        let Self {
            mut summary,
            items,
            zoom_items,
            ..
        } = self;

        debug_assert!(items.is_empty());
        for zoom_item in zoom_items.iter() {
            debug_assert!(zoom_item.live_info.is_none());
            debug_assert!(zoom_item.records.is_empty());
        }

        if summary.total_items == 0 {
            summary.min_val = 0.0;
            summary.max_val = 0.0;
        }
        ChromProcessedData(summary)
    }
}

impl ChromProcess for BigWigFullProcess {
    type Value = Value;
    async fn do_process<E: Error + Send + 'static>(
        &mut self,
        current_val: Value,
        next_val: Option<&Value>,
    ) -> Result<(), ProcessChromError<E>> {
        let Self {
            summary,
            items,
            zoom_items,
            ftx,
            chrom_id,
            options,
            runtime,
            chrom,
            length,
        } = self;
        let chrom_id = *chrom_id;
        let options = *options;
        let length = *length;

        BigWigWrite::process_val(
            current_val,
            next_val,
            length,
            &chrom,
            summary,
            items,
            options,
            &runtime,
            ftx,
            chrom_id,
        )
        .await?;

        BigWigWrite::process_val_zoom(
            zoom_items,
            options,
            current_val,
            next_val,
            &runtime,
            chrom_id,
        )
        .await?;

        Ok(())
    }
}

async fn encode_section(
    compress: bool,
    items_in_section: Vec<Value>,
    chrom_id: u32,
) -> io::Result<(SectionData, usize)> {
    use libdeflater::{CompressionLvl, Compressor};

    let mut bytes = Vec::with_capacity(24 + (items_in_section.len() * 24));

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
