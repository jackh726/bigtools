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

#[allow(clippy::too_many_arguments)]
pub fn intervals_to_array<R: BBIFileRead>(
    py: Python<'_>,
    bw: &mut BigWigRead<R>,
    chrom: &str,
    start: Option<i32>,
    end: Option<i32>,
    bins: Option<usize>,
    summary: Summary,
    exact: bool,
    uncovered: Option<f64>,
    oob: f64,
    fillna: Option<f64>,
    arr: Option<PyObject>,
) -> PyResult<PyObject> {
    let mut bbi = BBIRead::BigWig(bw);
    records_to_array(
        py, &mut bbi, chrom, start, end, bins, summary, exact, uncovered, oob, fillna, arr,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn entries_to_array<R: BBIFileRead>(
    py: Python<'_>,
    bb: &mut BigBedRead<R>,
    chrom: &str,
    start: Option<i32>,
    end: Option<i32>,
    bins: Option<usize>,
    summary: Summary,
    exact: bool,
    uncovered: Option<f64>,
    oob: f64,
    fillna: Option<f64>,
    arr: Option<PyObject>,
) -> PyResult<PyObject> {
    let mut bbi = BBIRead::BigBed(bb);
    records_to_array(
        py, &mut bbi, chrom, start, end, bins, summary, exact, uncovered, oob, fillna, arr,
    )
}

#[allow(clippy::too_many_arguments)]
fn records_to_array<R: BBIFileRead>(
    py: Python<'_>,
    bbi: &mut BBIRead<R>,
    chrom: &str,
    start: Option<i32>,
    end: Option<i32>,
    bins: Option<usize>,
    summary: Summary,
    exact: bool,
    uncovered: Option<f64>,
    oob: f64,
    fillna: Option<f64>,
    arr: Option<PyObject>,
) -> PyResult<PyObject> {
    let (start, end, chrom_length) = bbi.resolve_range(chrom, start, end)?;
    let (valid_start, valid_end) = (start.max(0) as u32, end.min(chrom_length) as u32);
    let arr = arr.unwrap_or_else(|| {
        let size = bins.unwrap_or((end - start) as usize);
        let init = uncovered.unwrap_or(f64::NAN);
        PyArray1::from_vec_bound(py, vec![init; size]).to_object(py)
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
                    fill_values(start, end, iter, uncovered, view).convert_err()?;
                }
                BBIRead::BigBed(bb) => {
                    let iter = bb
                        .get_interval(chrom, valid_start, valid_end)
                        .convert_err()?
                        .map(|item| item.map(BBIRecord::BedEntry));
                    fill_values(start, end, iter, uncovered, view).convert_err()?;
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
                        fill_binned(start, end, iter, summary, bins, uncovered, view)
                            .convert_err()?;
                    }
                    BBIRead::BigBed(bb) => {
                        let iter = bb
                            .get_zoom_interval(chrom, valid_start, valid_end, reduction_level)
                            .convert_err()?
                            .map(|item| item.map(BBIRecord::ZoomRecord));
                        fill_binned(start, end, iter, summary, bins, uncovered, view)
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
                        fill_binned(start, end, iter, summary, bins, uncovered, view)
                            .convert_err()?;
                    }
                    BBIRead::BigBed(bb) => {
                        let iter = bb
                            .get_interval(chrom, valid_start, valid_end)
                            .convert_err()?;
                        let coverage_iter = CoverageIterator::new(iter)
                            .convert_err()?
                            .map(|result| result.map(BBIRecord::Value));
                        fill_binned(start, end, coverage_iter, summary, bins, uncovered, view)
                            .convert_err()?;
                    }
                },
            }
        }
    };

    // Post-process fill for NaN positions. Applied before filling out-of-bounds.
    let mut array = array.readwrite();
    let mut view = array.as_array_mut();
    if let Some(fill) = fillna {
        for v in view.iter_mut() {
            if v.is_nan() {
                *v = fill;
            }
        }
    }

    // Handle out-of-bounds (bins/positions that overlap out-of-chromosome coordinates).
    fill_out_of_bounds(start, end, chrom_length, bins, oob, view.view_mut());

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
/// * `uncovered` - Value to use for bases with no data; `None` leaves them as NaN
/// * `view` - Mutable numpy array view to fill (must have length == end - start)
fn fill_values<I: Iterator<Item = Result<BBIRecord, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    uncovered: Option<f64>,
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

    // Finalize: replace remaining NaN with `uncovered` if a fill was requested.
    if let Some(fill) = uncovered {
        for val in view.iter_mut() {
            if val.is_nan() {
                *val = fill;
            }
        }
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
/// `uncovered` is the per-base value assigned to uncovered bases within each bin. `None`
/// means "exclude from the summary" and yields the classical statistics over only the
/// covered bases. `Some(v)` (commonly `Some(0.0)`) means "treat uncovered bases as having
/// value `v`"; passing `Some(0.0)` reproduces UCSC's `*0` summary statistics. The choice
/// applies uniformly: there are no separate `*0` variants of `Mean`, `Std`, `Min`, `Max`.
///
/// # Arguments
/// * `start` - Start coordinate of the target range
/// * `end` - End coordinate of the target range
/// * `iter` - Iterator of data intervals
/// * `summary` - How to aggregate values
/// * `bins` - Number of bins to divide the range into
/// * `uncovered` - Per-base value for uncovered bases; `None` excludes them
/// * `view` - Mutable numpy array view to fill (must have length == bins)
fn fill_binned<I: Iterator<Item = Result<BBIRecord, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    summary: Summary,
    bins: usize,
    uncovered: Option<f64>,
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
        for (i, stat) in stats
            .iter_mut()
            .enumerate()
            .take(last_bin + 1)
            .skip(first_bin)
        {
            let (bin_start, bin_end) = integer_bin_bounds(start, i, bin_size);

            // Calculate integer overlap between interval and bin
            let overlap_start = interval_start.max(bin_start);
            let overlap_end = interval_end.min(bin_end);
            let overlap = (overlap_end - overlap_start) as f64;

            if overlap > 0.0 {
                match interval {
                    BBIRecord::Value(v) => {
                        let value = v.value as f64;
                        stat.sum += value * overlap;
                        stat.sum_squares += value * value * overlap;
                        stat.min = stat.min.min(value);
                        stat.max = stat.max.max(value);
                        stat.bases_covered += overlap;
                    }
                    BBIRecord::ZoomRecord(z) => {
                        let summary = z.summary;
                        let overlap_factor = overlap / (z.end - z.start) as f64;
                        stat.sum += summary.sum * overlap_factor;
                        stat.sum_squares += summary.sum_squares * overlap_factor;
                        stat.min = stat.min.min(summary.min_val);
                        stat.max = stat.max.max(summary.max_val);
                        stat.bases_covered += summary.bases_covered as f64 * overlap_factor;
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

    // Fill the output array.
    //
    // `uncovered` is the per-base value assigned to uncovered bases within each bin:
    //   - `None`: uncovered bases are excluded from the summary computation
    //     ("regular" statistics over just the covered bases). For a fully-uncovered
    //     bin the result is NaN for stats over an empty set; for `BasesCovered` and
    //     `BinCovered` the result is the factual 0.
    //   - `Some(v)`: uncovered bases are folded into the summary as if they had
    //     value `v`. `Some(0.0)` reproduces UCSC's `*0` statistics (e.g.
    //     `mean0 = sum / bin_size`). Other values are well-defined generalizations.
    //
    // `BasesCovered` and `BinCovered` always report counts of *actual* data bases,
    // independent of `uncovered`.
    for i in 0..bins {
        let (bin_start, bin_end) = integer_bin_bounds(start, i, bin_size);
        let bin_width = (bin_end - bin_start) as f64;
        let n_covered = stats[i].bases_covered;
        let n_uncovered = bin_width - n_covered;

        // Effective sum/sum_squares and denominator after optionally folding in
        // `uncovered`. When `uncovered` is `None`, only covered bases contribute.
        let (s, ss, n) = match uncovered {
            None => (stats[i].sum, stats[i].sum_squares, n_covered),
            Some(u) => (
                stats[i].sum + u * n_uncovered,
                stats[i].sum_squares + u * u * n_uncovered,
                bin_width,
            ),
        };

        *view.index_mut(i) = match summary {
            Summary::Mean => {
                if n > 0.0 {
                    s / n
                } else {
                    f64::NAN
                }
            }
            Summary::Std => {
                if n < 1.0 {
                    f64::NAN
                } else {
                    let mut variance = ss - s * s / n;
                    if n > 1.0 {
                        variance /= n - 1.0;
                    }
                    variance.sqrt()
                }
            }
            Summary::Sum => {
                if n > 0.0 {
                    s
                } else {
                    f64::NAN
                }
            }
            Summary::SumSquares => {
                if n > 0.0 {
                    ss
                } else {
                    f64::NAN
                }
            }
            Summary::Min => match uncovered {
                None => {
                    if n_covered > 0.0 {
                        stats[i].min
                    } else {
                        f64::NAN
                    }
                }
                Some(u) if n_uncovered > 0.0 => {
                    if n_covered > 0.0 {
                        stats[i].min.min(u)
                    } else {
                        u
                    }
                }
                Some(_) => stats[i].min,
            },
            Summary::Max => match uncovered {
                None => {
                    if n_covered > 0.0 {
                        stats[i].max
                    } else {
                        f64::NAN
                    }
                }
                Some(u) if n_uncovered > 0.0 => {
                    if n_covered > 0.0 {
                        stats[i].max.max(u)
                    } else {
                        u
                    }
                }
                Some(_) => stats[i].max,
            },
            Summary::BasesCovered => n_covered,
            Summary::BinCovered => n_covered / bin_size,
        };
    }

    Ok(())
}

/// Compute the integer-aligned base coordinates of the bin at index `i`.
///
/// Float-valued `bin_size` does not generally fall on integer base boundaries, but the
/// underlying data is indexed by integer positions. We resolve this by truncating each
/// bin edge: bin `i` nominally spans
/// `[start + floor(i * bin_size), start + floor((i + 1) * bin_size))`.
///
/// Consequence: for non-integer `bin_size`, adjacent bins can have integer widths that
/// differ by 1 (e.g. for `bin_size = 2.4`, widths alternate 2, 2, 3, 2, 3) even though
/// each bin "intends" to span `bin_size` bases. The integer width returned here is the
/// actual number of base positions assigned to the bin.
///
/// For upsampling (`bin_size < 1`), the truncation can collapse a bin to zero integer
/// width. We expand such bins to span a single base so that overlap and division
/// operations stay well-defined; consecutive upsampled bins may then share the same
/// underlying base.
fn integer_bin_bounds(start: i32, i: usize, bin_size: f64) -> (i32, i32) {
    let bin_start = start + (i as f64 * bin_size) as i32;
    let bin_end = start + ((i + 1) as f64 * bin_size) as i32;
    if bin_start == bin_end {
        (bin_start, bin_start + 1)
    } else {
        (bin_start, bin_end)
    }
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
            None,
            arr.view_mut(),
        )
        .unwrap();

        // First bin (0-10): overlaps 5-10 = 5 bp with value 2.0
        // Second bin (10-20): overlaps 10-15 = 5 bp with value 2.0
        // uncovered=None excludes uncovered bases → mean of just the covered bp.
        assert_eq!(arr.to_vec(), vec![2.0, 2.0]);
    }

    #[test]
    fn test_fill_binned_partial_overlap() {
        // Bin 0 [0,10): overlaps [5,8) — 3 bp at value 4.0
        // Bin 1 [10,20): no overlap
        let start = 0;
        let end = 20;
        let n_bins = 2;

        // uncovered = None ("exclude uncovered" / classical statistics).
        // Bin 0: mean over the 3 covered bases = 4.0.
        // Bin 1: fully uncovered → NaN.
        let intervals = vec![Value {
            start: 5,
            end: 8,
            value: 4.0,
        }];
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            intervals.into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Mean,
            n_bins,
            None,
            arr.view_mut(),
        )
        .unwrap();
        assert!((arr[0] - 4.0).abs() < 1e-12);
        assert!(arr[1].is_nan());

        // uncovered = Some(0.0) (UCSC mean0 semantics: uncovered bases count as 0).
        // Bin 0: (4.0 * 3 + 0.0 * 7) / 10 = 1.2.
        // Bin 1: (0.0 * 10) / 10 = 0.0.
        let intervals = vec![Value {
            start: 5,
            end: 8,
            value: 4.0,
        }];
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            intervals.into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Mean,
            n_bins,
            Some(0.0),
            arr.view_mut(),
        )
        .unwrap();
        assert!((arr[0] - 1.2).abs() < 1e-12);
        assert!((arr[1] - 0.0).abs() < 1e-12);
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
            Some(-1.0),
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
            Some(-1.0),
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
            None,
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

    #[test]
    fn test_fill_binned_min_max_with_uncovered_zero_noninteger_binsize() {
        // 5 bins over range 0..12, bin_size = 2.4
        // Integer-aligned bin spans: [0,2), [2,4), [4,7), [7,9), [9,12)
        // Widths:                      2      2      3      2      3
        //
        // Single Value covering [0, 5) with value=5 covers bins 0 and 1 fully and
        // overlaps bin 2 partially (1 of 3 bases, at position 4). Bins 3 and 4 have
        // no coverage.
        //
        // With `uncovered = 0.0`, uncovered bases are folded into min/max as zeros.
        // The integer-adjusted bin width is load-bearing: comparing `bases_covered`
        // (=2 for bin 0) against the float `bin_size` (=2.4) would falsely report
        // bin 0 as partially covered, contaminating min/max with 0. Using the actual
        // integer bin width (=2) correctly identifies bin 0 as fully covered.
        let start = 0;
        let end = 12;
        let n_bins = 5;
        let uncovered = Some(0.0);

        // Positive value: min should reveal the "contamination" bug if present
        // (would pull min from 5 down to 0 in the fully-covered bins).
        let intervals = vec![Value {
            start: 0,
            end: 5,
            value: 5.0,
        }];
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            intervals.clone().into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Min,
            n_bins,
            uncovered,
            arr.view_mut(),
        )
        .unwrap();
        // Bin 0, 1 fully covered → min = 5. Bin 2 partial → min(5, 0) = 0.
        // Bins 3, 4 fully uncovered → uncovered (0).
        assert_eq!(arr.to_vec(), vec![5.0, 5.0, 0.0, 0.0, 0.0]);

        // Max with a negative value: partial coverage should pull max *up* to 0.
        let intervals_neg = vec![Value {
            start: 0,
            end: 5,
            value: -5.0,
        }];
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            intervals_neg.into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Max,
            n_bins,
            uncovered,
            arr.view_mut(),
        )
        .unwrap();
        // Bins 0, 1 fully covered → max = -5. Bin 2 partial → max(-5, 0) = 0.
        // Bins 3, 4 fully uncovered → uncovered (0).
        assert_eq!(arr.to_vec(), vec![-5.0, -5.0, 0.0, 0.0, 0.0]);

        // Max with positive value: partial coverage doesn't change max since
        // max(5, 0) == 5; the test pins that the integer-aligned check doesn't
        // accidentally pull max *down* in a fully-covered bin.
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            intervals.into_iter().map(BBIRecord::Value).map(Ok),
            Summary::Max,
            n_bins,
            uncovered,
            arr.view_mut(),
        )
        .unwrap();
        assert_eq!(arr.to_vec(), vec![5.0, 5.0, 5.0, 0.0, 0.0]);
    }

    /// Build an iterator of `Ok(BBIRecord::Value(...))` from a list of triples.
    fn vals(items: Vec<(u32, u32, f32)>) -> impl Iterator<Item = Result<BBIRecord, BBIReadError>> {
        items
            .into_iter()
            .map(|(s, e, v)| Value {
                start: s,
                end: e,
                value: v,
            })
            .map(BBIRecord::Value)
            .map(Ok)
    }

    #[test]
    fn test_fill_binned_sum_sum_squares() {
        // 2 bins over range 0..20. Values:
        //   [0, 5) value=2  → bin 0 contributes 2*5 = 10
        //   [5, 15) value=3 → bin 0 contributes 3*5 = 15; bin 1 contributes 3*5 = 15
        // Bin 0: sum=25 over 10 covered bp; bin 1: sum=15 over 5 covered bp.
        let start = 0;
        let end = 20;
        let n_bins = 2;
        let intervals = vec![(0u32, 5u32, 2.0f32), (5, 15, 3.0)];

        // Sum, uncovered=None: only covered contributes.
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            vals(intervals.clone()),
            Summary::Sum,
            n_bins,
            None,
            arr.view_mut(),
        )
        .unwrap();
        assert_eq!(arr.to_vec(), vec![25.0, 15.0]);

        // Sum, uncovered=Some(1.0): bin 1 also gets 1.0 * 5 uncovered bases.
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            vals(intervals.clone()),
            Summary::Sum,
            n_bins,
            Some(1.0),
            arr.view_mut(),
        )
        .unwrap();
        assert_eq!(arr.to_vec(), vec![25.0, 15.0 + 1.0 * 5.0]);

        // SumSquares, uncovered=None.
        // Bin 0: 2*2*5 + 3*3*5 = 20 + 45 = 65. Bin 1: 3*3*5 = 45.
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            vals(intervals.clone()),
            Summary::SumSquares,
            n_bins,
            None,
            arr.view_mut(),
        )
        .unwrap();
        assert_eq!(arr.to_vec(), vec![65.0, 45.0]);

        // SumSquares, uncovered=Some(2.0): adds 2*2 = 4 per uncovered base.
        // Bin 1 has 5 uncovered bp → +20.
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            vals(intervals),
            Summary::SumSquares,
            n_bins,
            Some(2.0),
            arr.view_mut(),
        )
        .unwrap();
        assert_eq!(arr.to_vec(), vec![65.0, 45.0 + 4.0 * 5.0]);
    }

    #[test]
    fn test_fill_binned_bases_and_bin_covered() {
        // 4 bins of width 5 over [0, 20). One value covering [0, 12).
        // Bin 0 [0, 5): fully covered (5).
        // Bin 1 [5, 10): fully covered (5).
        // Bin 2 [10, 15): partially covered (2 of 5).
        // Bin 3 [15, 20): no coverage (0).
        let start = 0;
        let end = 20;
        let n_bins = 4;
        let intervals = vec![(0u32, 12u32, 1.0f32)];

        // BasesCovered: actual integer count of covered bases per bin.
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            vals(intervals.clone()),
            Summary::BasesCovered,
            n_bins,
            None,
            arr.view_mut(),
        )
        .unwrap();
        assert_eq!(arr.to_vec(), vec![5.0, 5.0, 2.0, 0.0]);

        // `uncovered` must not affect BasesCovered.
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            vals(intervals.clone()),
            Summary::BasesCovered,
            n_bins,
            Some(99.0),
            arr.view_mut(),
        )
        .unwrap();
        assert_eq!(arr.to_vec(), vec![5.0, 5.0, 2.0, 0.0]);

        // BinCovered: bases_covered / bin_size.
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            vals(intervals),
            Summary::BinCovered,
            n_bins,
            None,
            arr.view_mut(),
        )
        .unwrap();
        assert_eq!(arr.to_vec(), vec![1.0, 1.0, 0.4, 0.0]);
    }

    #[test]
    fn test_fill_binned_std() {
        // 2 bins over [0, 10). Bin 0 sees two distinct values; bin 1 is empty.
        //   [0, 3) v=1  → contributes value 1 over 3 bases in bin 0
        //   [3, 5) v=4  → contributes value 4 over 2 bases in bin 0
        // Bin 0 sum = 1*3 + 4*2 = 11; sum_squares = 1*3 + 16*2 = 35; n_covered = 5.
        // Sample variance = (35 - 11*11/5) / (5-1) = (35 - 24.2) / 4 = 2.7
        // std = sqrt(2.7).
        let start = 0;
        let end = 10;
        let n_bins = 2;
        let intervals = vec![(0u32, 3u32, 1.0f32), (3, 5, 4.0)];

        // uncovered=None: classical sample std over only the covered bases.
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            vals(intervals.clone()),
            Summary::Std,
            n_bins,
            None,
            arr.view_mut(),
        )
        .unwrap();
        assert!((arr[0] - 2.7_f64.sqrt()).abs() < 1e-12);
        assert!(arr[1].is_nan()); // empty bin → NaN

        // uncovered=Some(0.0): the 0-padded computation.
        // Bin 0: s = 11 + 0*0 = 11; ss = 35 + 0 = 35; n = 5 (bin_width).
        //   variance = (35 - 121/5) / 4 = (35 - 24.2)/4 = 2.7
        //   std = sqrt(2.7). Same as None because bin 0 is fully covered.
        // Bin 1: s = 0; ss = 0; n = 5. variance = 0/4 = 0. std = 0.
        let mut arr = Array::from(vec![0.0; n_bins]);
        fill_binned(
            start,
            end,
            vals(intervals),
            Summary::Std,
            n_bins,
            Some(0.0),
            arr.view_mut(),
        )
        .unwrap();
        assert!((arr[0] - 2.7_f64.sqrt()).abs() < 1e-12);
        assert!((arr[1] - 0.0).abs() < 1e-12);
    }
}
