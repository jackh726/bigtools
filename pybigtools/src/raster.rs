use bigtools::{BBIFileRead, BBIReadError, BedEntry, BigBedRead, BigWigRead, Value, ZoomRecord};
use numpy::ndarray::ArrayViewMut;
use numpy::{PyArray1, PyArrayMethods};
use pyo3::exceptions;
use pyo3::prelude::*;
use std::ops::IndexMut;

use crate::coverage::CoverageIterator;
use crate::errors::ConvertResult;

#[derive(Copy, Clone, Debug)]
pub enum Summary {
    Mean,
    Std,
    Min,
    Max,
    Mean0,
    Std0,
    Min0,
    Max0,
    Sum,
    SumSquares,
    BasesCovered,
    BinCovered,
}

#[derive(Debug)]
enum BBIRecord {
    Value(Value),
    BedEntry(BedEntry),
    ZoomRecord(ZoomRecord),
}

#[derive(Copy, Clone, Debug)]
struct SummaryStats {
    sum: f64,
    sum_squares: f64,
    min: f64,
    max: f64,
    bases_covered: f64,
}

enum BBIRead<'a, R: BBIFileRead> {
    BigWig(&'a mut BigWigRead<R>),
    BigBed(&'a mut BigBedRead<R>),
}

impl<R: BBIFileRead> BBIRead<'_, R> {
    pub fn chroms(&self) -> &[bigtools::ChromInfo] {
        match self {
            BBIRead::BigWig(bw) => bw.chroms(),
            BBIRead::BigBed(bb) => bb.chroms(),
        }
    }

    /// Resolve a genomic range with potentially unspecified start and end positions.
    /// Returns a tuple of (start, end, chrom_length).
    pub fn resolve_range(
        &self,
        chrom_name: &str,
        start: Option<i32>,
        end: Option<i32>,
    ) -> PyResult<(i32, i32, i32)> {
        let chroms = self.chroms();
        let chrom = chroms.iter().find(|x| x.name == chrom_name);
        let length = match chrom {
            None => {
                return Err(PyErr::new::<exceptions::PyKeyError, _>(format!(
                    "No chromomsome with name `{}` found.",
                    chrom_name
                )))
            }
            Some(c) => c.length as i32,
        };
        Ok((start.unwrap_or(0), end.unwrap_or(length), length))
    }

    pub fn zoom_headers(&self) -> &[bigtools::ZoomHeader] {
        match self {
            BBIRead::BigWig(bw) => bw.info().zoom_headers.as_slice(),
            BBIRead::BigBed(bb) => bb.info().zoom_headers.as_slice(),
        }
    }

    /// Find the closest zoom resolution that does not exceed the desired bin size.
    /// If no such zoom level exists, returns None.
    pub fn best_zoom(&self, start: i32, end: i32, bins: usize) -> Option<u32> {
        let max_zoom_size = ((end - start) as f32 / (bins * 2) as f32) as u32;
        self.zoom_headers()
            .iter()
            .filter(|z| z.reduction_level <= max_zoom_size)
            .min_by_key(|z| max_zoom_size - z.reduction_level)
            .map(|z| z.reduction_level)
    }
}

pub fn intervals_to_array<R: BBIFileRead>(
    py: Python<'_>,
    bw: &mut BigWigRead<R>,
    chrom: &str,
    start: Option<i32>,
    end: Option<i32>,
    bins: Option<usize>,
    summary: Summary,
    exact: bool,
    missing: f64,
    oob: f64,
    arr: Option<PyObject>,
) -> PyResult<PyObject> {
    let mut bbi = BBIRead::BigWig(bw);
    records_to_array(
        py, &mut bbi, chrom, start, end, bins, summary, exact, missing, oob, arr,
    )
}

pub fn entries_to_array<R: BBIFileRead>(
    py: Python<'_>,
    bb: &mut BigBedRead<R>,
    chrom: &str,
    start: Option<i32>,
    end: Option<i32>,
    bins: Option<usize>,
    summary: Summary,
    exact: bool,
    missing: f64,
    oob: f64,
    arr: Option<PyObject>,
) -> PyResult<PyObject> {
    let mut bbi = BBIRead::BigBed(bb);
    records_to_array(
        py, &mut bbi, chrom, start, end, bins, summary, exact, missing, oob, arr,
    )
}

fn records_to_array<R: BBIFileRead>(
    py: Python<'_>,
    bbi: &mut BBIRead<R>,
    chrom: &str,
    start: Option<i32>,
    end: Option<i32>,
    bins: Option<usize>,
    summary: Summary,
    exact: bool,
    missing: f64,
    oob: f64,
    arr: Option<PyObject>,
) -> PyResult<PyObject> {
    let (start, end, chrom_length) = bbi.resolve_range(chrom, start, end)?;
    let (valid_start, valid_end) = (start.max(0) as u32, end.min(chrom_length) as u32);
    let arr = arr.unwrap_or_else(|| {
        let size = bins.unwrap_or((end - start) as usize);
        PyArray1::from_vec_bound(py, vec![missing; size]).to_object(py)
    });
    let array: &Bound<'_, PyArray1<f64>> =
        arr.downcast_bound::<PyArray1<f64>>(py).map_err(|_| {
            PyErr::new::<exceptions::PyValueError, _>(
                "`arr` option must be a one-dimensional numpy array, if passed.",
            )
        })?;

    // Fill the array
    match bins {
        // Base-level raster
        None => {
            if array.len()? != (end - start) as usize {
                return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
                    "Passed `arr` does not have the expected length (expected `{}`, found `{}`).",
                    (end - start) as usize,
                    array.len()?,
                )));
            }
            let mut array = array.readwrite();
            let view = array.as_array_mut();
            match bbi {
                BBIRead::BigWig(bw) => {
                    let iter = bw
                        .get_interval(chrom, valid_start, valid_end)
                        .convert_err()?
                        .map(|item| item.map(BBIRecord::Value));
                    fill_values(start, end, iter, missing, view).convert_err()?;
                }
                BBIRead::BigBed(bb) => {
                    let iter = bb
                        .get_interval(chrom, valid_start, valid_end)
                        .convert_err()?
                        .map(|item| item.map(BBIRecord::BedEntry));
                    fill_values(start, end, iter, missing, view).convert_err()?;
                }
            }
        }
        // Reduction
        Some(bins) => {
            if array.len()? != bins {
                return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
                    "Passed `arr` does not have the expected length (expected `{}`, found `{}`).",
                    bins,
                    array.len()?,
                )));
            }
            let reduction_level = if !exact {
                bbi.best_zoom(start, end, bins)
            } else {
                None
            };
            let mut array = array.readwrite();
            let view = array.as_array_mut();
            match reduction_level {
                // Interpolate from zoom-level summaries
                Some(reduction_level) => match bbi {
                    BBIRead::BigWig(bw) => {
                        let iter = bw
                            .get_zoom_interval(chrom, valid_start, valid_end, reduction_level)
                            .convert_err()?
                            .map(|item| item.map(BBIRecord::ZoomRecord));
                        fill_binned(start, end, iter, summary, bins, missing, view)
                            .convert_err()?;
                    }
                    BBIRead::BigBed(bb) => {
                        let iter = bb
                            .get_zoom_interval(chrom, valid_start, valid_end, reduction_level)
                            .convert_err()?
                            .map(|item| item.map(BBIRecord::ZoomRecord));
                        fill_binned(start, end, iter, summary, bins, missing, view)
                            .convert_err()?;
                    }
                },
                // Aggregate from base-level data
                None => match bbi {
                    BBIRead::BigWig(bw) => {
                        let iter = bw
                            .get_interval(chrom, valid_start, valid_end)
                            .convert_err()?
                            .map(|item| item.map(BBIRecord::Value));
                        fill_binned(start, end, iter, summary, bins, missing, view)
                            .convert_err()?;
                    }
                    BBIRead::BigBed(bb) => {
                        let iter = bb
                            .get_interval(chrom, valid_start, valid_end)
                            .convert_err()?;
                        let coverage_iter = CoverageIterator::new(iter)
                            .convert_err()?
                            .map(|result| result.map(BBIRecord::Value));
                        fill_binned(start, end, coverage_iter, summary, bins, missing, view)
                            .convert_err()?;
                    }
                },
            }
        }
    };

    // Handle out-of-bounds
    let mut array = array.readwrite();
    let view = array.as_array_mut();
    fill_out_of_bounds(start, end, chrom_length, bins, oob, view);

    Ok(arr)
}

/// Rasterize intervals at nucleotide resolution.
///
/// For [`Value`] intervals, the values are assigned over the covered bases.
/// For [`BedEntry`] intervals, the total coverage is accumulated over each base.
/// For [`ZoomRecord`] intervals, the bin mean statistic is assigned over the covered bases.
///
/// # Arguments
/// * `start` - Start coordinate of the target range
/// * `end` - End coordinate of the target range
/// * `iter` - Iterator of data intervals
/// * `missing` - Value to use for bases with no data
/// * `view` - Mutable numpy array view to fill (must have length == end - start)
fn fill_values<I: Iterator<Item = Result<BBIRecord, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    missing: f64,
    mut view: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(view.len(), (end - start) as usize);
    view.fill(f64::NAN);

    for interval in iter {
        let interval = interval?;

        // Clamp interval to our target range
        let interval_start = match interval {
            BBIRecord::Value(v) => (v.start as i32).max(start),
            BBIRecord::BedEntry(ref e) => (e.start as i32).max(start),
            BBIRecord::ZoomRecord(ref z) => (z.start as i32).max(start),
        };
        let interval_end = match interval {
            BBIRecord::Value(v) => (v.end as i32).min(end),
            BBIRecord::BedEntry(ref e) => (e.end as i32).min(end),
            BBIRecord::ZoomRecord(ref z) => (z.end as i32).min(end),
        };
        if interval_start >= interval_end {
            continue;
        }

        // Fill values
        let rel_start = (interval_start - start) as usize;
        let rel_end = (interval_end - start) as usize;
        let value = match interval {
            BBIRecord::Value(v) => v.value as f64,
            BBIRecord::BedEntry(_) => 1.0,
            BBIRecord::ZoomRecord(z) => z.summary.sum / (z.end - z.start) as f64,
        };
        for i in rel_start..rel_end {
            let val = *view.index_mut(i);
            *view.index_mut(i) = if val.is_nan() { value } else { val + value };
        }
    }

    // Finalize
    for val in view.iter_mut() {
        *val = if val.is_nan() { missing } else { *val };
    }

    Ok(())
}

/// Rasterize intervals over a uniform grid of bins.
///
/// This function iterates over each data interval and accumulates its contribution to the
/// overlapping bins, weighted by the overlap amount.
///
/// For [`Value`] intervals, the values are aggregated based on bin overlap.
/// For [`ZoomRecord`] intervals, the summary statistics are aggregated based on bin overlap.
///
/// # Notes
/// [`BedEntry`] intervals are not supported. They must be converted to coverage [`Value`]
/// intervals. See [`CoverageIterator`] for details.
///
/// Per-bin `bases_covered` is normalized up to the next integer (matching UCSC's convention),
/// with `sum` and `sum_squares` rescaled proportionally so that `mean = sum / bases_covered`
/// is invariant. This only has an effect on the [`ZoomRecord`] path, where partial overlap of
/// a zoom record can produce a fractional share of its `validCount`.
///
/// # Arguments
/// * `start` - Start coordinate of the target range
/// * `end` - End coordinate of the target range  
/// * `iter` - Iterator of data intervals
/// * `summary` - How to aggregate values
/// * `bins` - Number of bins to divide the range into
/// * `missing` - Value to use for bins with no intersecting data
/// * `view` - Mutable numpy array view to fill (must have length == bins)
fn fill_binned<I: Iterator<Item = Result<BBIRecord, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    summary: Summary,
    bins: usize,
    missing: f64,
    mut view: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(view.len(), bins);
    let bin_size = (end - start) as f64 / bins as f64;

    // Initialize statistics for each bin
    let mut stats = vec![
        SummaryStats {
            sum: 0.0,
            sum_squares: 0.0,
            min: f64::NAN,
            max: f64::NAN,
            bases_covered: 0.0,
        };
        bins
    ];

    // Process each data interval
    for interval_result in iter {
        let interval = interval_result?;

        // Clamp interval to our target range
        let interval_start = match interval {
            BBIRecord::Value(v) => (v.start as i32).max(start),
            BBIRecord::BedEntry(ref e) => (e.start as i32).max(start),
            BBIRecord::ZoomRecord(ref z) => (z.start as i32).max(start),
        };
        let interval_end = match interval {
            BBIRecord::Value(v) => (v.end as i32).min(end),
            BBIRecord::BedEntry(ref e) => (e.end as i32).min(end),
            BBIRecord::ZoomRecord(ref z) => (z.end as i32).min(end),
        };
        if interval_start >= interval_end {
            continue;
        }

        // Determine which bins this interval overlaps.
        // For the last bin, we need to handle the case where `interval_end` (open-ended) exactly
        // aligns with a bin boundary (in which case, we exclude the bin to the right).
        let rel_start = interval_start - start;
        let rel_end = interval_end - start;
        let first_bin = ((rel_start as f64) / bin_size).floor() as usize;
        let last_bin = if (rel_end as f64) % bin_size == 0.0 {
            (((rel_end as f64) / bin_size) - 1.0).floor() as usize
        } else {
            ((rel_end as f64) / bin_size).floor() as usize
        };
        // Clamp bins to `[0, bins - 1]` in case of floating point inaccuracies
        let first_bin = first_bin.min(bins - 1);
        let last_bin = last_bin.min(bins - 1);

        // Accumulate contribution of interval to each overlapping bin
        for i in first_bin..=last_bin {
            // Round bin edges down to the nearest integer (base) for overlap calculation
            let bin_start = start + (i as f64 * bin_size) as i32;
            let bin_end = start + ((i + 1) as f64 * bin_size) as i32;

            // For upsampling (case where bin_size < 1), ensure at least one base of overlap is counted
            let bin_end = if bin_start == bin_end {
                bin_start + 1
            } else {
                bin_end
            };

            // Calculate integer overlap between interval and adjusted bin
            let overlap_start = interval_start.max(bin_start);
            let overlap_end = interval_end.min(bin_end);
            let overlap = (overlap_end - overlap_start) as f64;

            if overlap > 0.0 {
                match interval {
                    BBIRecord::Value(v) => {
                        let value = v.value as f64;
                        stats[i].sum += value * overlap;
                        stats[i].sum_squares += value * value * overlap;
                        stats[i].min = stats[i].min.min(value);
                        stats[i].max = stats[i].max.max(value);
                        stats[i].bases_covered += overlap;
                    }
                    BBIRecord::ZoomRecord(z) => {
                        let summary = z.summary;
                        let overlap_factor = overlap / (z.end - z.start) as f64;
                        stats[i].sum += summary.sum * overlap_factor;
                        stats[i].sum_squares += summary.sum_squares * overlap_factor;
                        stats[i].min = stats[i].min.min(summary.min_val);
                        stats[i].max = stats[i].max.max(summary.max_val);
                        stats[i].bases_covered += summary.bases_covered as f64 * overlap_factor;
                    }
                    BBIRecord::BedEntry(_) => {
                        panic!("BED records must be converted to non-overlapping coverage intervals prior to binning.");
                    }
                }
            }
        }
    }

    // Normalize bases_covered to an integer count, matching UCSC's convention.
    //
    // For the Value path, accumulated `bases_covered` is already integer-valued (overlaps are
    // computed in integer coordinates), so `ceil` is a no-op.
    //
    // For the ZoomRecord path, a partial overlap of a zoom record contributes a fractional
    // share of its `validCount` (e.g. a zoom summarizing 10 bases with validCount=7, half
    // overlapped, contributes 3.5). UCSC rounds the per-bin total up to the next integer and
    // rescales `sum` and `sum_squares` by `ceil(bc) / bc` to keep `mean = sum / bases_covered`
    // invariant. The `ceil` direction preserves the invariant that `bases_covered >= 1`
    // whenever any data overlaps the bin. Note that `std` is *not* preserved under this
    // rescaling, but matching UCSC requires it.
    for stat in stats.iter_mut() {
        if stat.bases_covered > 0.0 {
            let vc = stat.bases_covered.ceil();
            let norm = vc / stat.bases_covered;
            stat.sum *= norm;
            stat.sum_squares *= norm;
            stat.bases_covered = vc;
        }
    }

    // Fill the output array
    for i in 0..bins {
        if stats[i].bases_covered == 0.0 {
            *view.index_mut(i) = missing;
        } else {
            match summary {
                Summary::Mean => {
                    *view.index_mut(i) = stats[i].sum / stats[i].bases_covered;
                }
                Summary::Std => {
                    let n = stats[i].bases_covered;
                    let mut variance = stats[i].sum_squares - (stats[i].sum * stats[i].sum) / n;
                    if n > 1.0 {
                        variance /= n - 1.0;
                    }
                    *view.index_mut(i) = variance.sqrt();
                }
                Summary::Min => {
                    *view.index_mut(i) = stats[i].min;
                }
                Summary::Max => {
                    *view.index_mut(i) = stats[i].max;
                }
                Summary::Mean0 => {
                    *view.index_mut(i) = stats[i].sum / bin_size;
                }
                Summary::Std0 => {
                    let n = bin_size;
                    let mut variance = stats[i].sum_squares - (stats[i].sum * stats[i].sum) / n;
                    if n > 1.0 {
                        variance /= n - 1.0;
                    }
                    *view.index_mut(i) = variance.sqrt();
                }
                Summary::Min0 => {
                    // Use integer-adjusted bin size when checking for partial coverage
                    let bin_start = start + (i as f64 * bin_size) as i32;
                    let bin_end = start + ((i + 1) as f64 * bin_size) as i32;
                    let adj_bin_size = (bin_end - bin_start) as f64;
                    let adj_bin_size = if adj_bin_size == 0.0 {
                        1.0
                    } else {
                        adj_bin_size
                    };
                    *view.index_mut(i) = if stats[i].bases_covered < adj_bin_size {
                        stats[i].min.min(0.0)
                    } else {
                        stats[i].min
                    }
                }
                Summary::Max0 => {
                    // Use integer-adjusted bin size when checking for partial coverage
                    let bin_start = start + (i as f64 * bin_size) as i32;
                    let bin_end = start + ((i + 1) as f64 * bin_size) as i32;
                    let adj_bin_size = (bin_end - bin_start) as f64;
                    let adj_bin_size = if adj_bin_size == 0.0 {
                        1.0
                    } else {
                        adj_bin_size
                    };
                    *view.index_mut(i) = if stats[i].bases_covered < adj_bin_size {
                        stats[i].max.max(0.0)
                    } else {
                        stats[i].max
                    }
                }
                Summary::Sum => {
                    *view.index_mut(i) = stats[i].sum;
                }
                Summary::SumSquares => {
                    *view.index_mut(i) = stats[i].sum_squares;
                }
                Summary::BasesCovered => {
                    *view.index_mut(i) = stats[i].bases_covered;
                }
                Summary::BinCovered => {
                    *view.index_mut(i) = stats[i].bases_covered / bin_size;
                }
            }
        }
    }

    Ok(())
}

/// Fill out-of-bounds elements in the array with a specified value.
///
/// If a bin overlaps any base coordinates outside the valid range of the chromosome, it is filled
/// with the `oob` value.
fn fill_out_of_bounds(
    start: i32,
    end: i32,
    chrom_length: i32,
    bins: Option<usize>,
    oob: f64,
    mut view: ArrayViewMut<'_, f64, numpy::Ix1>,
) {
    let bin_size = match bins {
        Some(bins) => (end - start) as f64 / bins as f64,
        _ => 1.0,
    };

    if start < 0 {
        for i in 0..view.len() {
            let bin_start = start + (i as f64 * bin_size) as i32;
            if bin_start < 0 {
                view[i] = oob;
            } else {
                break;
            }
        }
    }
    if end > chrom_length {
        for i in (0..view.len()).rev() {
            let bin_end = start + ((i + 1) as f64 * bin_size) as i32;
            if bin_end > chrom_length {
                view[i] = oob;
            } else {
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;
    use bigtools::Value;
    use numpy::ndarray::Array;

    #[test]
    fn test_fill_binned_single_interval() {
        let start = 0;
        let end = 20;
        let n_bins = 2;
        let intervals = vec![Value {
            start: 5,
            end: 15,
            value: 2.0,
        }];
        let mut arr = Array::from(vec![0.0; n_bins]);

        fill_binned(
            start,
            end,
            intervals.into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Mean,
            n_bins,
            0.0,
            arr.view_mut(),
        )
        .unwrap();

        // First bin (0-10): overlaps 5-10 = 5 bp with value 2.0
        // Second bin (10-20): overlaps 10-15 = 5 bp with value 2.0
        // Both should have value 2.0
        assert_eq!(arr.to_vec(), vec![2.0, 2.0]);
    }

    #[test]
    fn test_fill_binned_partial_overlap() {
        let start = 0;
        let end = 20;
        let n_bins = 2;
        let mut arr = Array::from(vec![0.0; n_bins]);

        let intervals = vec![Value {
            start: 5,
            end: 8,
            value: 4.0,
        }];
        fill_binned(
            start,
            end,
            intervals.into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Mean,
            n_bins,
            -1.0,
            arr.view_mut(),
        )
        .unwrap();

        // First bin (0-10): overlaps 5-8 = 3 bp with value 4.0
        // Second bin (10-20): no overlap
        assert_eq!(arr.to_vec(), vec![4.0, -1.0]);

        let intervals = vec![Value {
            start: 5,
            end: 8,
            value: 4.0,
        }];
        fill_binned(
            start,
            end,
            intervals.into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Mean0,
            n_bins,
            -1.0,
            arr.view_mut(),
        )
        .unwrap();

        // First bin (0-10): overlaps 5-8 = 3 bp with value 4.0
        // Second bin (10-20): no overlap
        assert!((arr[0] - (4.0 * 3.0 / 10.0)).abs() < 0.01);
        assert!((arr[1] - (-1.0)).abs() < 0.01);
    }

    #[test]
    fn test_fill_binned_min_max() {
        let start = 0;
        let end = 20;
        let n_bins = 2;
        let intervals = vec![
            Value {
                start: 0,
                end: 15,
                value: 1.0,
            },
            Value {
                start: 5,
                end: 25,
                value: 3.0,
            },
        ];

        // Test Min
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            intervals.clone().into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Min,
            n_bins,
            -1.0,
            arr.view_mut(),
        )
        .unwrap();

        // First bin (0-10): sees values 1.0 and 3.0, min = 1.0
        // Second bin (10-20): sees values 1.0 and 3.0, min = 1.0
        assert_eq!(arr.to_vec(), vec![1.0, 1.0]);

        // Test Max
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            intervals.into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Max,
            n_bins,
            -1.0,
            arr.view_mut(),
        )
        .unwrap();

        // First bin (0-10): sees values 1.0 and 3.0, max = 3.0
        // Second bin (10-20): sees values 1.0 and 3.0, max = 3.0
        assert_eq!(arr.to_vec(), vec![3.0, 3.0]);
    }

    #[test]
    fn test_fill_binned_noninteger_binsize() {
        // 5 bins over range 0..12
        // bin_size = 12 / 5 = 2.4
        // Bin edges are truncated to integers for overlap determination.
        //
        // Bin 0: [0, 2.4) → [0, 2)
        //      overlaps [1,3) by [1, 2) = 1 bp with value 1.0
        //      → 1.0
        // Bin 1: [2.4, 4.8) → [2, 4)
        //      overlaps [1,3) by [2, 3) = 1 bp with value 1.0
        //      → 1.0
        // Bin 2: [4.8, 7.2) → [4, 7)
        //      overlaps [4,8) by [4, 7) = 3 bp with value 2.0
        //      → 2.0
        // Bin 3: [7.2, 9.6) → [7, 9)
        //      overlaps [4,8) by [7, 8) = 1 bp with value 2.0
        //      overlaps [8,10) by [8, 9) = 1 bp with value 3.0
        //      → (2.0*1 + 3.0*1) / (1+1) = 2.5
        // Bin 4: [9.6, 12) → [9, 12)
        //      overlaps [8,10) by [9, 10) = 1 bp with value 3.0
        //      → 3.0
        let start = 0;
        let end = 12;
        let n_bins = 5;
        let intervals = vec![
            Value {
                start: 1,
                end: 3,
                value: 1.0,
            },
            Value {
                start: 4,
                end: 8,
                value: 2.0,
            },
            Value {
                start: 8,
                end: 10,
                value: 3.0,
            },
        ];
        let answer = vec![1.0, 1.0, 2.0, 2.5, 3.0];
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            intervals.into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Mean,
            n_bins,
            0.0,
            arr.view_mut(),
        )
        .unwrap();

        for i in 0..n_bins {
            assert!(
                (arr[i] - answer[i]).abs() < 0.01,
                "Bin {}: expected {}, got {}, diff {}",
                i,
                answer[i],
                arr[i],
                (arr[i] - answer[i]).abs()
            );
        }
    }
}
