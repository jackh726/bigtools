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
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::future::Future;
use std::io::{self, BufWriter, Write};

use futures::future::FutureExt;
use futures::sink::SinkExt;

use byteorder::{NativeEndian, WriteBytesExt};
use tokio::runtime::{Handle, Runtime};

use crate::utils::chromvalues::ChromValues;
use crate::utils::tell::Tell;
use crate::{write_info, ChromData, ChromProcessingInputSectionChannel};

use crate::bbi::{Summary, Value, ZoomRecord, BIGWIG_MAGIC};
use crate::bbiwrite::{
    self, encode_zoom_section, write_blank_headers, write_zooms, BBIWriteOptions,
    ProcessChromError, SectionData,
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

    /// Write the values from `V` as a bigWig. Will utilize the provided runtime for encoding values and for reading through the values (potentially parallelized by chromosome).
    pub fn write<
        Values: ChromValues<Value = Value> + Send + 'static,
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
            let fut = BigWigWrite::process_chrom(
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
        Values: ChromValues<Value = Value>,
        V: ChromData<Values = Values>,
    >(
        self,
        chrom_sizes: HashMap<String, u32>,
        vals: V,
        runtime: Runtime,
    ) -> Result<(), ProcessChromError<Values::Error>> {
        self.write_internal(chrom_sizes, vals, runtime, BigWigWrite::process_chrom)
    }

    fn write_internal<
        Values: ChromValues<Value = Value>,
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

        let (total_summary_offset, full_data_offset, pre_data) = BigWigWrite::write_pre(&mut file)?;

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
        let (data_size, chrom_index_start, index_start, total_sections) = bbiwrite::write_mid(
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
            uncompress_buf_size,
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
        let (data_size, chrom_index_start, index_start, total_sections) = bbiwrite::write_mid(
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

    async fn process_val<I: ChromValues<Value = Value>>(
        current_val: Value,
        chrom_length: u32,
        chrom: &String,
        chrom_values: &mut I,
        summary: &mut Summary,
        items: &mut Vec<Value>,
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
        match chrom_values.peek() {
            None | Some(Err(_)) => (),
            Some(Ok(next_val)) => {
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
        if chrom_values.peek().is_none() || items.len() >= options.items_per_slot as usize {
            let items = std::mem::take(items);
            let handle = runtime
                .spawn(encode_section(options.compress, items, chrom_id))
                .map(|f| f.unwrap());
            ftx.send(handle.boxed()).await.expect("Couldn't send");
        }

        Ok(())
    }

    async fn process_val_zoom<I: ChromValues<Value = Value>>(
        zoom_items: &mut Vec<ZoomItem>,
        options: BBIWriteOptions,
        current_val: Value,
        chrom_values: &mut I,
        runtime: &Handle,
        chrom_id: u32,
    ) -> Result<(), ProcessChromError<I::Error>> {
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
                    && chrom_values.peek().is_none()
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
                    if chrom_values.peek().is_none() {
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

    pub(crate) async fn process_chrom<I: ChromValues<Value = Value>>(
        zooms_channels: Vec<(u32, ChromProcessingInputSectionChannel)>,
        mut ftx: ChromProcessingInputSectionChannel,
        chrom_id: u32,
        options: BBIWriteOptions,
        runtime: Handle,
        mut chrom_values: I,
        chrom: String,
        chrom_length: u32,
    ) -> Result<Summary, ProcessChromError<I::Error>> {
        let mut summary = Summary {
            total_items: 0,
            bases_covered: 0,
            min_val: f64::MAX,
            max_val: f64::MIN,
            sum: 0.0,
            sum_squares: 0.0,
        };

        let mut items = Vec::with_capacity(options.items_per_slot as usize);
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

            BigWigWrite::process_val(
                current_val,
                chrom_length,
                &chrom,
                &mut chrom_values,
                &mut summary,
                &mut items,
                options,
                &runtime,
                &mut ftx,
                chrom_id,
            )
            .await?;

            BigWigWrite::process_val_zoom(
                &mut zoom_items,
                options,
                current_val,
                &mut chrom_values,
                &runtime,
                chrom_id,
            )
            .await?;
        }

        debug_assert!(items.is_empty());
        for zoom_item in zoom_items.iter_mut() {
            debug_assert!(zoom_item.live_info.is_none());
            debug_assert!(zoom_item.records.is_empty());
        }

        if summary.total_items == 0 {
            summary.min_val = 0.0;
            summary.max_val = 0.0;
        }
        Ok(summary)
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

            BigWigWrite::process_val(
                current_val,
                chrom_length,
                &chrom,
                &mut chrom_values,
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

            BigWigWrite::process_val_zoom(
                &mut zoom_items,
                options,
                current_val,
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
