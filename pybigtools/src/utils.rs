use std::collections::VecDeque;
use std::ops::IndexMut;

#[cfg(feature = "remote")]
use bigtools::utils::file::remote_file::RemoteFile;
use bigtools::utils::file::reopen::ReopenableFile;
use bigtools::{
    BBIFileRead, BBIReadError, BedEntry, BigBedRead, BigWigRead, CachedBBIFileRead, Value,
    ZoomRecord,
};
use numpy::ndarray::ArrayViewMut;
use numpy::{PyArray1, PyArrayMethods};
use pyo3::exceptions;
use pyo3::prelude::*;

use crate::errors::ConvertResult;
use crate::file_like::PyFileLikeObject;

pub enum BBIReadRaw {
    Closed,
    BigWigFile(BigWigRead<CachedBBIFileRead<ReopenableFile>>),
    #[cfg(feature = "remote")]
    BigWigRemote(BigWigRead<CachedBBIFileRead<RemoteFile>>),
    BigWigFileLike(BigWigRead<CachedBBIFileRead<PyFileLikeObject>>),
    BigBedFile(BigBedRead<CachedBBIFileRead<ReopenableFile>>),
    #[cfg(feature = "remote")]
    BigBedRemote(BigBedRead<CachedBBIFileRead<RemoteFile>>),
    BigBedFileLike(BigBedRead<CachedBBIFileRead<PyFileLikeObject>>),
}

#[derive(Copy, Clone, Debug)]
pub enum Summary {
    Mean,
    Min,
    Max,
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
    let (start, end, chrom_length) = bigwig_start_end_length(bw, chrom, start, end)?;
    let (intervals_start, intervals_end) = (start.max(0) as u32, end.min(chrom_length) as u32);
    let arr = arr.unwrap_or_else(|| {
        let size = bins.unwrap_or((end - start) as usize);
        PyArray1::from_vec_bound(py, vec![missing; size]).to_object(py)
    });

    // Fill array
    let array: &Bound<'_, PyArray1<f64>> =
        arr.downcast_bound::<PyArray1<f64>>(py).map_err(|_| {
            PyErr::new::<exceptions::PyValueError, _>(
                "`arr` option must be a one-dimensional numpy array, if passed.",
            )
        })?;
    match bins {
        // Base-level data
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
            let iter = bw
                .get_interval(&chrom, intervals_start, intervals_end)
                .convert_err()?;
            fill_values(start, end, iter, missing, view).convert_err()?;
        }
        // Binned data
        Some(bins) => {
            if array.len()? != bins {
                return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
                    "Passed `arr` does not have the expected length (expected `{}`, found `{}`).",
                    bins,
                    array.len()?,
                )));
            }
            let zoom = if !exact {
                bigwig_closest_zoom(bw, start, end, bins)
            } else {
                None
            };
            let mut array = array.readwrite();
            let view = array.as_array_mut();
            match zoom {
                // Interpolate from zoom-level summaries
                Some(zoom) => {
                    let iter = bw
                        .get_zoom_interval(
                            &chrom,
                            intervals_start,
                            intervals_end,
                            zoom.reduction_level,
                        )
                        .convert_err()?;
                    fill_binned_values_from_zoom(start, end, iter, summary, bins, missing, view)
                        .convert_err()?;
                }
                // Aggregate from base-level data
                None => {
                    let iter = bw
                        .get_interval(&chrom, intervals_start, intervals_end)
                        .convert_err()?;
                    fill_binned_values_from_intervals(
                        start, end, iter, summary, bins, missing, view,
                    )
                    .convert_err()?;
                }
            }
        }
    };

    // Handle out-of-bounds
    let mut array = array.readwrite();
    let view = array.as_array_mut();
    fill_out_of_bounds(start, end, chrom_length, bins, oob, view);

    Ok(arr)
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
    let (start, end, length) = bigbed_start_end_length(bb, chrom, start, end)?;
    let (intervals_start, intervals_end) = (start.max(0) as u32, end.min(length) as u32);
    let arr = arr.unwrap_or_else(|| {
        let size = bins.unwrap_or((end - start) as usize);
        PyArray1::from_vec_bound(py, vec![missing; size]).to_object(py)
    });

    // Fill array
    let array: &Bound<'_, PyArray1<f64>> =
        arr.downcast_bound::<PyArray1<f64>>(py).map_err(|_| {
            PyErr::new::<exceptions::PyValueError, _>(
                "Passed `arr` option must be a one-dimensional numpy array.",
            )
        })?;
    match bins {
        // Base-level data
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
            let iter = bb
                .get_interval(&chrom, intervals_start, intervals_end)
                .convert_err()?;
            fill_coverage(start, end, iter, missing, view).convert_err()?;
        }
        // Binned data
        Some(bins) => {
            if array.len()? != bins {
                return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
                    "Passed `arr` does not have the expected length (expected `{}`, found `{}`).",
                    bins,
                    array.len()?,
                )));
            }
            let zoom = if !exact {
                bigbed_closest_zoom(bb, start, end, bins)
            } else {
                None
            };
            let mut array = array.readwrite();
            let view = array.as_array_mut();
            match zoom {
                // Interpolate from zoom-level summaries
                Some(zoom) => {
                    let iter = bb
                        .get_zoom_interval(
                            &chrom,
                            intervals_start,
                            intervals_end,
                            zoom.reduction_level,
                        )
                        .convert_err()?;
                    fill_binned_coverage_from_zoom(start, end, iter, summary, bins, missing, view)
                        .convert_err()?;
                }
                // Aggregate from base-level data
                None => {
                    let iter = bb
                        .get_interval(&chrom, intervals_start, intervals_end)
                        .convert_err()?;
                    fill_binned_coverage_from_entries(
                        start, end, iter, summary, bins, missing, view,
                    )
                    .convert_err()?;
                }
            }
        }
    };

    // Handle out-of-bounds
    let mut array = array.readwrite();
    let view = array.as_array_mut();
    fill_out_of_bounds(start, end, length, bins, oob, view);

    Ok(arr)
}

fn fill_values<I: Iterator<Item = Result<Value, BBIReadError>>>(
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
        let interval_start = ((interval.start as i32) - start) as usize;
        let interval_end = ((interval.end as i32) - start) as usize;
        for i in interval_start..interval_end {
            let val = *view.index_mut(i);
            *view.index_mut(i) = if val.is_nan() {
                interval.value as f64
            } else {
                val + interval.value as f64
            };
        }
    }
    for val in view.iter_mut() {
        *val = if val.is_nan() { missing } else { *val };
    }
    Ok(())
}

fn fill_coverage<I: Iterator<Item = Result<BedEntry, BBIReadError>>>(
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
        let interval_start = ((interval.start as i32) - start) as usize;
        let interval_end = ((interval.end as i32) - start) as usize;
        for i in interval_start..interval_end {
            let val = *view.index_mut(i);
            *view.index_mut(i) = if val.is_nan() { 1.0 } else { val + 1.0 };
        }
    }
    for val in view.iter_mut() {
        *val = if val.is_nan() { missing } else { *val };
    }

    Ok(())
}

fn fill_binned_values_from_intervals<I: Iterator<Item = Result<Value, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    summary: Summary,
    bins: usize,
    missing: f64,
    mut view: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(view.len(), bins);
    view.fill(missing);

    // (bin, bin_start, bin_end, Option<(covered_bases, value)>)
    let mut deque: VecDeque<(usize, i32, i32, Option<(i32, f64)>)> = VecDeque::new();
    let bin_size = (end - start) as f64 / bins as f64;
    for interval in iter {
        let interval = interval?;
        let interval_start = (interval.start as i32).max(start) - start;
        let interval_end = (interval.end as i32).min(end) - start;
        let bin_start = ((interval_start as f64) / bin_size) as usize;
        let bin_end = (((interval_end - 1) as f64) / bin_size) as usize;

        while let Some(front) = deque.front_mut() {
            if front.0 < bin_start {
                let front = deque.pop_front().unwrap();
                let bin = front.0;

                match summary {
                    Summary::Min => {
                        view[bin] = front.3.map(|v| v.1).unwrap_or(missing);
                    }
                    Summary::Max => {
                        view[bin] = front.3.map(|v| v.1).unwrap_or(missing);
                    }
                    Summary::Mean => {
                        view[bin] = front.3.map(|(c, v)| v / (c as f64)).unwrap_or(missing);
                    }
                }
            } else {
                break;
            }
        }
        while let Some(bin) = deque
            .back()
            .map(|last_bin| {
                if last_bin.0 < bin_end {
                    Some(last_bin.0 + 1)
                } else {
                    None
                }
            })
            .unwrap_or(Some(bin_start))
        {
            let bin_start = ((bin as f64) * bin_size) as i32;
            let bin_end = (((bin + 1) as f64) * bin_size) as i32;

            deque.push_back((bin, bin_start, bin_end, None));
        }
        for bin in bin_start..bin_end {
            assert!(deque.iter().find(|b| b.0 == bin).is_some());
        }
        for (_bin, bin_start, bin_end, data) in deque.iter_mut() {
            // Since bed files can have overlapping entries, a previous entry could have extended
            // past the last bin for the current interval
            // Not relevant here, but keeping for parity
            if interval_end <= *bin_start {
                break;
            }
            let (c, v) = data.get_or_insert_with(|| {
                match summary {
                    // min & max are defined for NAN and we are about to set it
                    // can't use 0.0 because it may be either below or above the real value
                    Summary::Min | Summary::Max => (0, f64::NAN),
                    // addition is not defined for NAN
                    Summary::Mean => (0, 0.0),
                }
            });
            match summary {
                Summary::Min => {
                    *v = v.min(interval.value as f64);
                }
                Summary::Max => {
                    *v = v.max(interval.value as f64);
                }
                Summary::Mean => {
                    let overlap_start = (*bin_start).max(interval_start);
                    let overlap_end = (*bin_end).min(interval_end);
                    let overlap_size: i32 = overlap_end - overlap_start;
                    *v += (overlap_size as f64) * interval.value as f64;
                    *c += overlap_size;
                }
            }
        }
    }
    while let Some(front) = deque.pop_front() {
        let bin = front.0;

        match summary {
            Summary::Min => {
                view[bin] = front.3.map(|v| v.1).unwrap_or(missing);
            }
            Summary::Max => {
                view[bin] = front.3.map(|v| v.1).unwrap_or(missing);
            }
            Summary::Mean => {
                view[bin] = front.3.map(|(c, v)| v / (c as f64)).unwrap_or(missing);
            }
        }
    }
    Ok(())
}

fn fill_binned_coverage_from_entries<I: Iterator<Item = Result<BedEntry, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    summary: Summary,
    bins: usize,
    missing: f64,
    mut view: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(view.len(), bins);
    view.fill(missing);

    // (<bin>, <bin_start>, <bin_end>, <covered_bases>, <sum>)
    // covered_bases = 0 if uncovered, 1 if covered
    let mut deque: VecDeque<(usize, i32, i32, Vec<i32>, Vec<f64>)> = VecDeque::new();
    let bin_size = (end - start) as f64 / bins as f64;
    for interval in iter {
        let interval = interval?;
        let interval_start = (interval.start as i32).max(start) - start;
        let interval_end = (interval.end as i32).min(end) - start;
        let bin_start = ((interval_start as f64) / bin_size) as usize;
        let bin_end = (((interval_end - 1) as f64) / bin_size) as usize;

        while let Some(front) = deque.front_mut() {
            if front.0 < bin_start {
                let front = deque.pop_front().unwrap();
                let bin = front.0;

                match summary {
                    Summary::Min => {
                        view[bin] = front
                            .4
                            .into_iter()
                            .reduce(|min, v| min.min(v))
                            .unwrap_or(missing);
                    }
                    Summary::Max => {
                        view[bin] = front
                            .4
                            .into_iter()
                            .reduce(|max, v| max.max(v))
                            .unwrap_or(missing);
                    }
                    Summary::Mean => {
                        view[bin] = front
                            .3
                            .iter()
                            .any(|v| *v > 0)
                            .then(|| front.3.into_iter().sum::<i32>() as f64)
                            // Map nan to 0 so the sum is non-nan
                            .map(|c| front.4.into_iter().map(|c| c.max(0.0)).sum::<f64>() / c)
                            .unwrap_or(missing);
                    }
                }
            } else {
                break;
            }
        }
        while let Some(bin) = deque
            .back()
            .map(|last_bin| {
                if last_bin.0 < bin_end {
                    Some(last_bin.0 + 1)
                } else {
                    None
                }
            })
            .unwrap_or(Some(bin_start))
        {
            let bin_start = ((bin as f64) * bin_size) as i32;
            let bin_end = (((bin + 1) as f64) * bin_size) as i32;

            deque.push_back((
                bin,
                bin_start,
                bin_end,
                vec![0; (bin_end - bin_start) as usize],
                vec![missing; (bin_end - bin_start) as usize],
            ));
        }
        for bin in bin_start..bin_end {
            assert!(deque.iter().find(|b| b.0 == bin).is_some());
        }
        for (_bin, bin_start, bin_end, covered, data) in deque.iter_mut() {
            // Since bed files can have overlapping entries, a previous entry could have extended
            // past the last bin for the current interval
            if interval_end <= *bin_start {
                break;
            }
            let overlap_start = (*bin_start).max(interval_start);
            let overlap_end = (*bin_end).min(interval_end);

            let range =
                ((overlap_start - *bin_start) as usize)..((overlap_end - *bin_start) as usize);
            for i in &mut data[range.clone()] {
                // If NAN, then 0.0 + 1.0, else i + 1.0
                *i = (*i).max(0.0) + 1.0;
            }
            for i in &mut covered[range] {
                *i = (*i).max(1);
            }
        }
    }
    while let Some(front) = deque.pop_front() {
        let bin = front.0;

        match summary {
            Summary::Min => {
                view[bin] = front
                    .4
                    .into_iter()
                    .reduce(|min, v| min.min(v))
                    .unwrap_or(missing);
            }
            Summary::Max => {
                view[bin] = front
                    .4
                    .into_iter()
                    .reduce(|max, v| max.max(v))
                    .unwrap_or(missing);
            }
            Summary::Mean => {
                view[bin] = front
                    .3
                    .iter()
                    .any(|v| *v > 0)
                    .then(|| front.3.into_iter().sum::<i32>() as f64)
                    // Map nan to 0 so the sum is non-nan
                    .map(|c| front.4.into_iter().map(|c| c.max(0.0)).sum::<f64>() / c)
                    .unwrap_or(missing);
            }
        }
    }
    Ok(())
}

fn fill_binned_values_from_zoom<I: Iterator<Item = Result<ZoomRecord, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    summary: Summary,
    bins: usize,
    missing: f64,
    mut view: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(view.len(), bins);
    view.fill(missing);

    // (bin, bin_start, bin_end, Option<(covered_bases, value)>)
    let mut deque: VecDeque<(usize, i32, i32, Option<(i32, f64)>)> = VecDeque::new();
    let bin_size = (end - start) as f64 / bins as f64;
    for interval in iter {
        let interval = interval?;
        let interval_start = (interval.start as i32).max(start) - start;
        let interval_end = (interval.end as i32).min(end) - start;
        let bin_start = ((interval_start as f64) / bin_size) as usize;
        let bin_end = (((interval_end - 1) as f64) / bin_size) as usize;

        while let Some(front) = deque.front_mut() {
            if front.0 < bin_start {
                let front = deque.pop_front().unwrap();
                let bin = front.0;

                match summary {
                    Summary::Min => {
                        view[bin] = front.3.map(|v| v.1).unwrap_or(missing);
                    }
                    Summary::Max => {
                        view[bin] = front.3.map(|v| v.1).unwrap_or(missing);
                    }
                    Summary::Mean => {
                        view[bin] = front.3.map(|(c, v)| v / (c as f64)).unwrap_or(missing);
                    }
                }
            } else {
                break;
            }
        }
        while let Some(bin) = deque
            .back()
            .map(|last_bin| {
                if last_bin.0 < bin_end {
                    Some(last_bin.0 + 1)
                } else {
                    None
                }
            })
            .unwrap_or(Some(bin_start))
        {
            let bin_start = ((bin as f64) * bin_size) as i32;
            let bin_end = (((bin + 1) as f64) * bin_size) as i32;

            deque.push_back((bin, bin_start, bin_end, None));
        }
        for bin in bin_start..bin_end {
            assert!(deque.iter().find(|b| b.0 == bin).is_some());
        }
        for (_bin, bin_start, bin_end, data) in deque.iter_mut() {
            // Since bed files can have overlapping entries, a previous entry could have extended
            // past the last bin for the current interval
            // Not relevant here, but keeping for parity
            if interval_end <= *bin_start {
                break;
            }
            let overlap_start = (*bin_start).max(interval_start);
            let overlap_end = (*bin_end).min(interval_end);
            let overlap_size: i32 = overlap_end - overlap_start;
            match data.as_mut() {
                Some((c, v)) => match summary {
                    Summary::Min => {
                        *v = v.min(interval.summary.min_val as f64);
                        *c += overlap_size;
                    }
                    Summary::Max => {
                        *v = v.max(interval.summary.max_val as f64);
                        *c += overlap_size;
                    }
                    Summary::Mean => {
                        let zoom_mean =
                            (interval.summary.sum as f64) / (interval.summary.bases_covered as f64);
                        *v += (overlap_size as f64) * zoom_mean;
                        *c += overlap_size;
                    }
                },
                None => match summary {
                    Summary::Min => {
                        *data = Some((overlap_size, interval.summary.min_val as f64));
                    }
                    Summary::Max => {
                        *data = Some((overlap_size, interval.summary.max_val as f64));
                    }
                    Summary::Mean => {
                        let zoom_mean =
                            (interval.summary.sum as f64) / (interval.summary.bases_covered as f64);
                        *data = Some((overlap_size, (overlap_size as f64) * zoom_mean));
                    }
                },
            }
        }
    }
    while let Some(front) = deque.pop_front() {
        let bin = front.0;

        match summary {
            Summary::Min => {
                view[bin] = front.3.map(|v| v.1).unwrap_or(missing);
            }
            Summary::Max => {
                view[bin] = front.3.map(|v| v.1).unwrap_or(missing);
            }
            Summary::Mean => {
                view[bin] = front.3.map(|(c, v)| v / (c as f64)).unwrap_or(missing);
            }
        }
    }
    Ok(())
}

fn fill_binned_coverage_from_zoom<I: Iterator<Item = Result<ZoomRecord, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    summary: Summary,
    bins: usize,
    missing: f64,
    mut view: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(view.len(), bins);
    view.fill(missing);

    // (<bin>, <bin_start>, <bin_end>, <covered_bases>, <sum>)
    // covered_bases = 0 if uncovered, 1 if covered
    let mut deque: VecDeque<(usize, i32, i32, Vec<i32>, Vec<f64>)> = VecDeque::new();
    let bin_size = (end - start) as f64 / bins as f64;
    for interval in iter {
        let interval = interval?;
        let interval_start = (interval.start as i32).max(start) - start;
        let interval_end = (interval.end as i32).min(end) - start;
        let bin_start = ((interval_start as f64) / bin_size) as usize;
        let bin_end = (((interval_end - 1) as f64) / bin_size) as usize;

        while let Some(front) = deque.front_mut() {
            if front.0 < bin_start {
                let front = deque.pop_front().unwrap();
                let bin = front.0;

                match summary {
                    Summary::Min => {
                        view[bin] = front
                            .4
                            .into_iter()
                            .reduce(|min, v| min.min(v))
                            .unwrap_or(missing);
                    }
                    Summary::Max => {
                        view[bin] = front
                            .4
                            .into_iter()
                            .reduce(|max, v| max.max(v))
                            .unwrap_or(missing);
                    }
                    Summary::Mean => {
                        view[bin] = front
                            .3
                            .iter()
                            .any(|v| *v > 0)
                            .then(|| front.3.into_iter().sum::<i32>() as f64)
                            // Map nan to 0 so the sum is non-nan
                            .map(|c| front.4.into_iter().map(|c| c.max(0.0)).sum::<f64>() / c)
                            .unwrap_or(missing);
                    }
                }
            } else {
                break;
            }
        }
        while let Some(bin) = deque
            .back()
            .map(|last_bin| {
                if last_bin.0 < bin_end {
                    Some(last_bin.0 + 1)
                } else {
                    None
                }
            })
            .unwrap_or(Some(bin_start))
        {
            let bin_start = ((bin as f64) * bin_size) as i32;
            let bin_end = (((bin + 1) as f64) * bin_size) as i32;

            deque.push_back((
                bin,
                bin_start,
                bin_end,
                vec![0; (bin_end - bin_start) as usize],
                vec![missing; (bin_end - bin_start) as usize],
            ));
        }
        for bin in bin_start..bin_end {
            assert!(deque.iter().find(|b| b.0 == bin).is_some());
        }
        for (_bin, bin_start, bin_end, covered, data) in deque.iter_mut() {
            // Since bed files can have overlapping entries, a previous entry could have extended
            // past the last bin for the current interval
            if interval_end <= *bin_start {
                break;
            }
            let overlap_start = (*bin_start).max(interval_start);
            let overlap_end = (*bin_end).min(interval_end);

            let mean = interval.summary.sum / (interval.summary.bases_covered) as f64;
            let range =
                ((overlap_start - *bin_start) as usize)..((overlap_end - *bin_start) as usize);
            for i in &mut data[range.clone()] {
                *i = i.max(0.0);
                match summary {
                    Summary::Mean => *i += mean,
                    Summary::Min => *i = i.min(interval.summary.min_val),
                    Summary::Max => *i = i.max(interval.summary.max_val),
                }
            }

            for i in &mut covered[range] {
                *i = (*i).max(1);
            }
        }
    }
    while let Some(front) = deque.pop_front() {
        let bin = front.0;

        match summary {
            Summary::Min => {
                view[bin] = front
                    .4
                    .into_iter()
                    .reduce(|min, v| min.min(v))
                    .unwrap_or(missing);
            }
            Summary::Max => {
                view[bin] = front
                    .4
                    .into_iter()
                    .reduce(|max, v| max.max(v))
                    .unwrap_or(missing);
            }
            Summary::Mean => {
                view[bin] = front
                    .3
                    .iter()
                    .any(|v| *v > 0)
                    .then(|| front.3.into_iter().sum::<i32>() as f64)
                    // Map nan to 0 so the sum is non-nan
                    .map(|c| front.4.into_iter().map(|c| c.max(0.0)).sum::<f64>() / c)
                    .unwrap_or(missing);
            }
        }
    }
    Ok(())
}

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
        let bin_start = 0;
        let interval_end = 0 - start;
        let bin_end = ((interval_end as f64) / bin_size).ceil() as usize;
        for i in bin_start..bin_end {
            view[i] = oob;
        }
    }
    if end > chrom_length {
        let interval_start = (chrom_length as i32) - start;
        let interval_end = end - start;
        let bin_start = ((interval_start as f64) / bin_size) as usize;
        let bin_end = ((interval_end as f64) / bin_size).ceil() as usize;
        assert_eq!(bin_end, view.len());
        for i in bin_start..bin_end {
            view[i] = oob;
        }
    }
}

fn bigwig_start_end_length<R>(
    bw: &BigWigRead<R>,
    chrom_name: &str,
    start: Option<i32>,
    end: Option<i32>,
) -> PyResult<(i32, i32, i32)> {
    let chroms = bw.chroms();
    start_end_length_inner(chrom_name, chroms, start, end)
}

fn bigbed_start_end_length<R>(
    bb: &BigBedRead<R>,
    chrom_name: &str,
    start: Option<i32>,
    end: Option<i32>,
) -> PyResult<(i32, i32, i32)> {
    let chroms = bb.chroms();
    start_end_length_inner(chrom_name, chroms, start, end)
}

fn start_end_length_inner(
    chrom_name: &str,
    chroms: &[bigtools::ChromInfo],
    start: Option<i32>,
    end: Option<i32>,
) -> PyResult<(i32, i32, i32)> {
    let chrom = chroms.into_iter().find(|x| x.name == chrom_name);
    let length = match chrom {
        None => {
            return Err(PyErr::new::<exceptions::PyKeyError, _>(format!(
                "No chromomsome with name `{}` found.",
                chrom_name
            )))
        }
        Some(c) => c.length as i32,
    };
    return Ok((start.unwrap_or(0), end.unwrap_or(length), length));
}

fn bigwig_closest_zoom<'a, R: BBIFileRead>(
    bw: &'a BigWigRead<R>,
    start: i32,
    end: i32,
    bins: usize,
) -> Option<&'a bigtools::ZoomHeader> {
    let max_zoom_size = ((end - start) as f32 / (bins * 2) as f32) as u32;
    bw.info()
        .zoom_headers
        .iter()
        .filter(|z| z.reduction_level <= max_zoom_size)
        .min_by_key(|z| max_zoom_size - z.reduction_level)
}

fn bigbed_closest_zoom<'a, R: BBIFileRead>(
    bb: &'a BigBedRead<R>,
    start: i32,
    end: i32,
    bins: usize,
) -> Option<&'a bigtools::ZoomHeader> {
    let max_zoom_size = ((end - start) as f32 / (bins * 2) as f32) as u32;
    bb.info()
        .zoom_headers
        .iter()
        .filter(|z| z.reduction_level <= max_zoom_size)
        .min_by_key(|z| max_zoom_size - z.reduction_level)
}

#[cfg(test)]
mod test {
    use bigtools::BedEntry;
    use numpy::ndarray::Array;

    use crate::utils::fill_binned_coverage_from_entries;
    use crate::utils::Summary;

    #[test]
    fn test_fill_binned_coverage_from_entries() {
        fn eq_with_nan_eq(a: f64, b: f64) -> bool {
            (a.is_nan() && b.is_nan()) || (a == b)
        }

        #[track_caller]
        fn assert_equal(found: Vec<f64>, expected: Vec<f64>) {
            let equal = (found.len() == expected.len()) &&  // zip stops at the shortest
            found.iter()
               .zip(expected.iter())
               .all(|(a,b)| eq_with_nan_eq(*a, *b));
            if !equal {
                panic!("Vecs not equal: expected {:?}, found {:?}", expected, found);
            }
        }

        let entries = [];
        let mut arr = Array::from(vec![f64::NAN]);
        fill_binned_coverage_from_entries(
            10,
            20,
            entries.into_iter().map(|v| Ok(v)),
            Summary::Mean,
            1,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [f64::NAN].into_iter().collect::<Vec<_>>());

        let entries = [];
        let mut arr = Array::from(vec![f64::NAN]);
        fill_binned_coverage_from_entries(
            10,
            20,
            entries.into_iter().map(|v| Ok(v)),
            Summary::Mean,
            1,
            0.0,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [0.0].into_iter().collect::<Vec<_>>());

        let entries = [BedEntry {
            start: 10,
            end: 20,
            rest: "".to_string(),
        }];
        let mut arr = Array::from(vec![f64::NAN]);
        fill_binned_coverage_from_entries(
            10,
            20,
            entries.into_iter().map(|v| Ok(v)),
            Summary::Mean,
            1,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [1.0].into_iter().collect::<Vec<_>>());

        let entries = [BedEntry {
            start: 10,
            end: 20,
            rest: "".to_string(),
        }];
        let mut arr = Array::from(vec![f64::NAN, f64::NAN]);
        fill_binned_coverage_from_entries(
            10,
            20,
            entries.into_iter().map(|v| Ok(v)),
            Summary::Mean,
            2,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [1.0, 1.0].into_iter().collect::<Vec<_>>());

        let entries = [
            BedEntry {
                start: 10,
                end: 20,
                rest: "".to_string(),
            },
            BedEntry {
                start: 15,
                end: 20,
                rest: "".to_string(),
            },
        ];
        let mut arr = Array::from(vec![f64::NAN, f64::NAN]);
        fill_binned_coverage_from_entries(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            Summary::Mean,
            2,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [1.0, 2.0].into_iter().collect::<Vec<_>>());
        let mut arr = Array::from(vec![f64::NAN]);
        fill_binned_coverage_from_entries(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            Summary::Min,
            1,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [1.0].into_iter().collect::<Vec<_>>());
        let mut arr = Array::from(vec![f64::NAN]);
        fill_binned_coverage_from_entries(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            Summary::Max,
            1,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [2.0].into_iter().collect::<Vec<_>>());

        let entries = [
            BedEntry {
                start: 10,
                end: 20,
                rest: "".to_string(),
            },
            BedEntry {
                start: 15,
                end: 20,
                rest: "".to_string(),
            },
            BedEntry {
                start: 15,
                end: 25,
                rest: "".to_string(),
            },
        ];
        let mut arr = Array::from(vec![f64::NAN, f64::NAN]);
        fill_binned_coverage_from_entries(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            Summary::Mean,
            2,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [1.0, 3.0].into_iter().collect::<Vec<_>>());
        let mut arr = Array::from(vec![f64::NAN]);
        fill_binned_coverage_from_entries(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            Summary::Min,
            1,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [1.0].into_iter().collect::<Vec<_>>());
        let mut arr = Array::from(vec![f64::NAN]);
        fill_binned_coverage_from_entries(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            Summary::Max,
            1,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [3.0].into_iter().collect::<Vec<_>>());
    }
}
