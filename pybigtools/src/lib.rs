#![allow(non_snake_case)]

use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufReader};
use std::ops::IndexMut;
use std::path::Path;

use bigtools::bed::autosql::parse::parse_autosql;
use bigtools::bed::bedparser::BedParser;
use bigtools::bedchromdata::BedParserStreamingIterator;
#[cfg(feature = "remote")]
use bigtools::utils::file::remote_file::RemoteFile;
use bigtools::utils::file::reopen::ReopenableFile;
use bigtools::utils::misc::{
    bigwig_average_over_bed, BigWigAverageOverBedEntry, BigWigAverageOverBedError, Name,
};
use bigtools::{
    BBIFileRead, BBIReadError, BedEntry, BigBedRead as BigBedReadRaw,
    BigBedWrite as BigBedWriteRaw, BigWigRead as BigWigReadRaw, BigWigWrite as BigWigWriteRaw,
    CachedBBIFileRead, Value, ZoomRecord,
};

use bigtools::utils::reopen::Reopen;
use file_like::PyFileLikeObject;
use numpy::ndarray::ArrayViewMut;
use numpy::PyArray1;
use pyo3::exceptions::{self, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyAny, PyDict, PyFloat, PyInt, PyIterator, PyString, PyTuple};
use pyo3::wrap_pyfunction;

use tokio::runtime;
use url::Url;

mod file_like;

type ValueTuple = (u32, u32, f32);

fn start_end(
    bbi: &BBIReadRaw,
    chrom_name: &str,
    start: Option<u32>,
    end: Option<u32>,
) -> PyResult<(u32, u32)> {
    let chroms = match &bbi {
        BBIReadRaw::BigWigFile(b) => b.chroms(),
        #[cfg(feature = "remote")]
        BBIReadRaw::BigWigRemote(b) => b.chroms(),
        BBIReadRaw::BigWigFileLike(b) => b.chroms(),
        BBIReadRaw::BigBedFile(b) => b.chroms(),
        #[cfg(feature = "remote")]
        BBIReadRaw::BigBedRemote(b) => b.chroms(),
        BBIReadRaw::BigBedFileLike(b) => b.chroms(),
    };
    let chrom = chroms.into_iter().find(|x| x.name == chrom_name);
    let length = match chrom {
        None => {
            return Err(PyErr::new::<exceptions::PyException, _>(format!(
                "No chromomsome with name `{}` found.",
                chrom_name
            )))
        }
        Some(c) => c.length,
    };
    return Ok((start.unwrap_or(0), end.unwrap_or(length)));
}

fn bigwig_start_end_length<R>(
    bbi: &BigWigReadRaw<R>,
    chrom_name: &str,
    start: Option<i32>,
    end: Option<i32>,
) -> PyResult<(i32, i32, i32)> {
    let chroms = bbi.chroms();
    start_end_length_inner(chrom_name, chroms, start, end)
}

fn bigbed_start_end_length<R>(
    bbi: &BigBedReadRaw<R>,
    chrom_name: &str,
    start: Option<i32>,
    end: Option<i32>,
) -> PyResult<(i32, i32, i32)> {
    let chroms = bbi.chroms();
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
            return Err(PyErr::new::<exceptions::PyException, _>(format!(
                "No chromomsome with name `{}` found.",
                chrom_name
            )))
        }
        Some(c) => c.length as i32,
    };
    return Ok((start.unwrap_or(0), end.unwrap_or(length), length));
}

impl Reopen for PyFileLikeObject {
    fn reopen(&self) -> io::Result<Self> {
        Ok(self.clone())
    }
}

enum BBIReadRaw {
    BigWigFile(BigWigReadRaw<CachedBBIFileRead<ReopenableFile>>),
    #[cfg(feature = "remote")]
    BigWigRemote(BigWigReadRaw<CachedBBIFileRead<RemoteFile>>),
    BigWigFileLike(BigWigReadRaw<CachedBBIFileRead<PyFileLikeObject>>),
    BigBedFile(BigBedReadRaw<CachedBBIFileRead<ReopenableFile>>),
    #[cfg(feature = "remote")]
    BigBedRemote(BigBedReadRaw<CachedBBIFileRead<RemoteFile>>),
    BigBedFileLike(BigBedReadRaw<CachedBBIFileRead<PyFileLikeObject>>),
}

#[derive(Copy, Clone, Debug)]
enum Summary {
    Mean,
    Min,
    Max,
}

trait ToPyErr {
    fn to_py_err(self) -> PyErr;
}

impl ToPyErr for bigtools::ZoomIntervalError {
    fn to_py_err(self) -> PyErr {
        PyErr::new::<exceptions::PyException, _>(format!("{}", self))
    }
}
impl ToPyErr for bigtools::BBIReadError {
    fn to_py_err(self) -> PyErr {
        PyErr::new::<exceptions::PyException, _>(format!("{}", self))
    }
}

trait ConvertResult<T> {
    fn convert_err(self) -> Result<T, PyErr>;
}
impl<T, E: ToPyErr> ConvertResult<T> for Result<T, E> {
    fn convert_err(self) -> Result<T, PyErr> {
        self.map_err(|e| e.to_py_err())
    }
}

fn intervals_to_array<R: BBIFileRead>(
    py: Python<'_>,
    b: &mut BigWigReadRaw<R>,
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
    let (start, end, length) = bigwig_start_end_length(b, chrom, start, end)?;
    let zoom = if let (Some(bins), false) = (bins, exact) {
        let max_zoom_size = ((end - start) as f32 / (bins * 2) as f32) as u32;
        let zoom = b
            .info()
            .zoom_headers
            .iter()
            .filter(|z| z.reduction_level <= max_zoom_size)
            .min_by_key(|z| max_zoom_size - z.reduction_level);
        zoom
    } else {
        None
    };
    let arr = match (bins, arr) {
        (_, Some(arr)) => arr,
        (Some(bins), None) => PyArray1::from_vec(py, vec![missing; bins]).to_object(py),
        (None, None) => PyArray1::from_vec(py, vec![missing; (end - start) as usize]).to_object(py),
    };
    let v: &PyArray1<f64> = arr.downcast::<PyArray1<f64>>(py).map_err(|_| {
        PyErr::new::<exceptions::PyException, _>(
            "`arr` option must be a one-dimensional numpy array, if passed.",
        )
    })?;
    let (intervals_start, intervals_end) = (start.max(0) as u32, end.min(length) as u32);
    let bin_size = match bins {
        Some(bins) => {
            if v.len() != bins {
                return Err(PyErr::new::<exceptions::PyException, _>(format!(
                    "`arr` does not the expected size (expected `{}`, found `{}`), if passed.",
                    bins,
                    v.len(),
                )));
            }

            let mut array = v.readwrite();

            match zoom {
                Some(zoom) => {
                    let iter = b
                        .get_zoom_interval(
                            &chrom,
                            intervals_start,
                            intervals_end,
                            zoom.reduction_level,
                        )
                        .convert_err()?;
                    to_array_zoom(
                        start,
                        end,
                        iter,
                        summary,
                        bins,
                        missing,
                        array.as_array_mut(),
                    )
                    .convert_err()?;
                }
                None => {
                    let iter = b
                        .get_interval(&chrom, intervals_start, intervals_end)
                        .convert_err()?;
                    to_array_bins(
                        start,
                        end,
                        iter,
                        summary,
                        bins,
                        missing,
                        array.as_array_mut(),
                    )
                    .convert_err()?;
                }
            };

            (end - start) as f64 / bins as f64
        }
        _ => {
            if v.len() != (end - start) as usize {
                return Err(PyErr::new::<exceptions::PyException, _>(format!(
                    "`arr` does not the expected size (expected `{}`, found `{}`), if passed.",
                    (end - start) as usize,
                    v.len(),
                )));
            }

            let mut array = v.readwrite();

            let iter = b
                .get_interval(&chrom, intervals_start, intervals_end)
                .convert_err()?;
            to_array(start, end, iter, missing, array.as_array_mut()).convert_err()?;

            (end - start) as f64
        }
    };

    {
        let mut array = v.readwrite();
        let mut array: ArrayViewMut<'_, f64, numpy::Ix1> = array.as_array_mut();

        if start < 0 {
            let bin_start = 0;
            let interval_end = 0 - start;
            let bin_end = ((interval_end as f64) / bin_size).ceil() as usize;
            for i in bin_start..bin_end {
                array[i] = oob;
            }
        }
        if end >= length {
            let interval_start = (length as i32) - start;
            let interval_end = end - start;
            let bin_start = ((interval_start as f64) / bin_size) as usize;
            let bin_end = ((interval_end as f64) / bin_size).ceil() as usize;
            assert_eq!(bin_end, array.len() - 1);
            for i in bin_start..bin_end {
                array[i] = oob;
            }
        }
    }

    Ok(arr)
}
fn entries_to_array<R: BBIFileRead>(
    py: Python<'_>,
    b: &mut BigBedReadRaw<R>,
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
    let (start, end, length) = bigbed_start_end_length(b, chrom, start, end)?;
    let zoom = if let (Some(bins), false) = (bins, exact) {
        let max_zoom_size = ((end - start) as f32 / (bins * 2) as f32) as u32;
        let zoom = b
            .info()
            .zoom_headers
            .iter()
            .filter(|z| z.reduction_level <= max_zoom_size)
            .min_by_key(|z| max_zoom_size - z.reduction_level);
        zoom
    } else {
        None
    };
    let arr = match (bins, arr) {
        (_, Some(arr)) => arr,
        (Some(bins), None) => PyArray1::from_vec(py, vec![missing; bins]).to_object(py),
        (None, None) => PyArray1::from_vec(py, vec![missing; (end - start) as usize]).to_object(py),
    };
    let v: &PyArray1<f64> = arr.downcast::<PyArray1<f64>>(py).map_err(|_| {
        PyErr::new::<exceptions::PyException, _>(
            "`arr` option must be a one-dimensional numpy array, if passed.",
        )
    })?;
    let (intervals_start, intervals_end) = (start.max(0) as u32, end.min(length) as u32);
    let bin_size = match bins {
        Some(bins) => {
            if v.len() != bins {
                return Err(PyErr::new::<exceptions::PyException, _>(format!(
                    "`arr` does not the expected size (expected `{}`, found `{}`), if passed.",
                    bins,
                    v.len(),
                )));
            }

            let mut array = v.readwrite();

            match zoom {
                Some(zoom) => {
                    let iter = b
                        .get_zoom_interval(
                            &chrom,
                            intervals_start,
                            intervals_end,
                            zoom.reduction_level,
                        )
                        .convert_err()?;
                    to_entry_array_zoom(
                        start,
                        end,
                        iter,
                        summary,
                        bins,
                        missing,
                        array.as_array_mut(),
                    )
                    .convert_err()?;
                }
                None => {
                    let iter = b
                        .get_interval(&chrom, intervals_start, intervals_end)
                        .convert_err()?;
                    to_entry_array_bins(
                        start,
                        end,
                        iter,
                        summary,
                        bins,
                        missing,
                        array.as_array_mut(),
                    )
                    .convert_err()?;
                }
            };

            (end - start) as f64 / bins as f64
        }
        _ => {
            if v.len() != (end - start) as usize {
                return Err(PyErr::new::<exceptions::PyException, _>(format!(
                    "`arr` does not the expected size (expected `{}`, found `{}`), if passed.",
                    (end - start) as usize,
                    v.len(),
                )));
            }

            let mut array = v.readwrite();

            let iter = b
                .get_interval(&chrom, intervals_start, intervals_end)
                .convert_err()?;
            to_entry_array(start, end, iter, missing, array.as_array_mut()).convert_err()?;

            (end - start) as f64
        }
    };

    {
        let mut array = v.readwrite();
        let mut array: ArrayViewMut<'_, f64, numpy::Ix1> = array.as_array_mut();

        if start < 0 {
            let bin_start = 0;
            let interval_end = 0 - start;
            let bin_end = ((interval_end as f64) / bin_size).ceil() as usize;
            for i in bin_start..bin_end {
                array[i] = oob;
            }
        }
        if end >= length {
            let interval_start = (length as i32) - start;
            let interval_end = end - start;
            let bin_start = ((interval_start as f64) / bin_size) as usize;
            let bin_end = ((interval_end as f64) / bin_size).ceil() as usize;
            assert_eq!(bin_end, array.len() - 1);
            for i in bin_start..bin_end {
                array[i] = oob;
            }
        }
    }

    Ok(arr)
}

fn to_array<I: Iterator<Item = Result<Value, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    missing: f64,
    mut v: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(v.len(), (end - start) as usize);
    v.fill(missing);
    for interval in iter {
        let interval = interval?;
        let interval_start = ((interval.start as i32) - start) as usize;
        let interval_end = ((interval.end as i32) - start) as usize;
        for i in interval_start..interval_end {
            *v.index_mut(i) += interval.value as f64;
        }
    }
    Ok(())
}
fn to_array_bins<I: Iterator<Item = Result<Value, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    summary: Summary,
    bins: usize,
    missing: f64,
    mut v: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(v.len(), bins);
    v.fill(missing);

    let mut bin_data: VecDeque<(usize, i32, i32, Option<f64>)> = VecDeque::new();
    let bin_size = (end - start) as f64 / bins as f64;
    for interval in iter {
        let interval = interval?;
        let interval_start = (interval.start as i32).max(start) - start;
        let interval_end = (interval.end as i32).min(end) - start;
        let bin_start = ((interval_start as f64) / bin_size) as usize;
        let bin_end = ((interval_end as f64) / bin_size).ceil() as usize;

        while let Some(front) = bin_data.front_mut() {
            if front.0 < bin_start {
                let front = bin_data.pop_front().unwrap();
                let bin = front.0;

                match summary {
                    Summary::Min => {
                        v[bin] = front.3.unwrap_or(missing);
                    }
                    Summary::Max => {
                        v[bin] = front.3.unwrap_or(missing);
                    }
                    Summary::Mean => {
                        v[bin] = front.3.map(|v| v / bin_size).unwrap_or(missing);
                    }
                }
            } else {
                break;
            }
        }
        while let Some(bin) = bin_data
            .back_mut()
            .map(|b| ((b.0 + 1) < bin_end).then(|| b.0 + 1))
            .unwrap_or(Some(bin_start))
        {
            let bin_start = ((bin as f64) * bin_size).ceil() as i32;
            let bin_end = (((bin + 1) as f64) * bin_size).ceil() as i32;

            bin_data.push_back((bin, bin_start, bin_end, None));
        }
        for bin in bin_data.iter() {
            assert!(
                (bin_start..bin_end).contains(&bin.0),
                "{} not in {}..{}",
                bin.0,
                bin_start,
                bin_end
            );
        }
        for bin in bin_start..bin_end {
            assert!(bin_data.iter().find(|b| b.0 == bin).is_some());
        }
        for (_bin, bin_start, bin_end, data) in bin_data.iter_mut() {
            let v = data.get_or_insert_with(|| {
                match summary {
                    // min & max are defined for NAN and we are about to set it
                    // can't use 0.0 because it may be either below or above the real value
                    Summary::Min | Summary::Max => f64::NAN,
                    // addition is not defined for NAN
                    Summary::Mean => 0.0,
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
                }
            }
        }
    }
    while let Some(front) = bin_data.pop_front() {
        let bin = front.0;

        match summary {
            Summary::Min => {
                v[bin] = front.3.unwrap_or(missing);
            }
            Summary::Max => {
                v[bin] = front.3.unwrap_or(missing);
            }
            Summary::Mean => {
                v[bin] = front.3.map(|v| v / bin_size).unwrap_or(missing);
            }
        }
    }
    Ok(())
}
fn to_array_zoom<I: Iterator<Item = Result<ZoomRecord, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    summary: Summary,
    bins: usize,
    missing: f64,
    mut v: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(v.len(), bins);
    v.fill(missing);

    let mut bin_data: VecDeque<(usize, i32, i32, Option<f64>)> = VecDeque::new();
    let bin_size = (end - start) as f64 / bins as f64;
    for interval in iter {
        let interval = interval?;
        let interval_start = (interval.start as i32).max(start) - start;
        let interval_end = (interval.end as i32).min(end) - start;
        let bin_start = ((interval_start as f64) / bin_size) as usize;
        let bin_end = ((interval_end as f64) / bin_size).ceil() as usize;

        while let Some(front) = bin_data.front_mut() {
            if front.0 < bin_start {
                let front = bin_data.pop_front().unwrap();
                let bin = front.0;

                match summary {
                    Summary::Min => {
                        v[bin] = front.3.unwrap_or(missing);
                    }
                    Summary::Max => {
                        v[bin] = front.3.unwrap_or(missing);
                    }
                    Summary::Mean => {
                        v[bin] = front.3.map(|v| v / bin_size).unwrap_or(missing);
                    }
                }
            } else {
                break;
            }
        }
        while let Some(bin) = bin_data
            .back_mut()
            .map(|b| ((b.0 + 1) < bin_end).then(|| b.0 + 1))
            .unwrap_or(Some(bin_start))
        {
            let bin_start = ((bin as f64) * bin_size).ceil() as i32;
            let bin_end = (((bin + 1) as f64) * bin_size).ceil() as i32;

            bin_data.push_back((bin, bin_start, bin_end, None));
        }
        for bin in bin_data.iter() {
            assert!(
                (bin_start..bin_end).contains(&bin.0),
                "{} not in {}..{}",
                bin.0,
                bin_start,
                bin_end
            );
        }
        for bin in bin_start..bin_end {
            assert!(bin_data.iter().find(|b| b.0 == bin).is_some());
        }
        for (_bin, bin_start, bin_end, data) in bin_data.iter_mut() {
            match data.as_mut() {
                Some(v) => match summary {
                    Summary::Min => {
                        *v = v.min(interval.summary.min_val as f64);
                    }
                    Summary::Max => {
                        *v = v.max(interval.summary.max_val as f64);
                    }
                    Summary::Mean => {
                        let overlap_start = (*bin_start).max(interval_start);
                        let overlap_end = (*bin_end).min(interval_end);
                        let overlap_size = overlap_end - overlap_start;
                        let zoom_mean =
                            (interval.summary.sum as f64) / (interval.summary.bases_covered as f64);
                        *v += (overlap_size as f64) * zoom_mean;
                    }
                },
                None => match summary {
                    Summary::Min => {
                        *data = Some(interval.summary.min_val as f64);
                    }
                    Summary::Max => {
                        *data = Some(interval.summary.max_val as f64);
                    }
                    Summary::Mean => {
                        let overlap_start = (*bin_start).max(interval_start);
                        let overlap_end = (*bin_end).min(interval_end);
                        let overlap_size = overlap_end - overlap_start;
                        let zoom_mean =
                            (interval.summary.sum as f64) / (interval.summary.bases_covered as f64);
                        *data = Some((overlap_size as f64) * zoom_mean);
                    }
                },
            }
        }
    }
    while let Some(front) = bin_data.pop_front() {
        let bin = front.0;

        match summary {
            Summary::Min => {
                v[bin] = front.3.unwrap_or(missing);
            }
            Summary::Max => {
                v[bin] = front.3.unwrap_or(missing);
            }
            Summary::Mean => {
                v[bin] = front.3.map(|v| v / bin_size).unwrap_or(missing);
            }
        }
    }
    Ok(())
}
fn to_entry_array<I: Iterator<Item = Result<BedEntry, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    missing: f64,
    mut v: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(v.len(), (end - start) as usize);
    v.fill(missing);
    for interval in iter {
        let interval = interval?;
        let interval_start = ((interval.start as i32) - start) as usize;
        let interval_end = ((interval.end as i32) - start) as usize;
        for i in interval_start..interval_end {
            *v.index_mut(i) += 1.0;
        }
    }
    Ok(())
}
fn to_entry_array_bins<I: Iterator<Item = Result<BedEntry, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    summary: Summary,
    bins: usize,
    missing: f64,
    mut v: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(v.len(), bins);
    v.fill(missing);

    let mut bin_data: VecDeque<(usize, i32, i32, Vec<f64>)> = VecDeque::new();
    let bin_size = (end - start) as f64 / bins as f64;
    for interval in iter {
        let interval = interval?;
        let interval_start = (interval.start as i32).max(start) - start;
        let interval_end = (interval.end as i32).min(end) - start;
        let bin_start = ((interval_start as f64) / bin_size) as usize;
        let bin_end = ((interval_end as f64) / bin_size).ceil() as usize;

        while let Some(front) = bin_data.front_mut() {
            if front.0 < bin_start {
                let front = bin_data.pop_front().unwrap();
                let bin = front.0;

                match summary {
                    Summary::Min => {
                        v[bin] = front
                            .3
                            .into_iter()
                            .reduce(|min, v| min.min(v))
                            .unwrap_or(missing);
                    }
                    Summary::Max => {
                        v[bin] = front
                            .3
                            .into_iter()
                            .reduce(|max, v| max.max(v))
                            .unwrap_or(missing);
                    }
                    Summary::Mean => {
                        v[bin] = front.3.into_iter().sum::<f64>() / bin_size;
                    }
                }
            } else {
                break;
            }
        }
        while let Some(bin) = bin_data
            .back_mut()
            .map(|b| ((b.0 + 1) < bin_end).then(|| b.0 + 1))
            .unwrap_or(Some(bin_start))
        {
            let bin_start = ((bin as f64) * bin_size).ceil() as i32;
            let bin_end = (((bin + 1) as f64) * bin_size).ceil() as i32;

            bin_data.push_back((
                bin,
                bin_start,
                bin_end,
                vec![missing; (bin_end - bin_start) as usize],
            ));
        }
        for bin in bin_data.iter() {
            assert!(
                (bin_start..bin_end).contains(&bin.0),
                "{} not in {}..{}",
                bin.0,
                bin_start,
                bin_end
            );
        }
        for bin in bin_start..bin_end {
            assert!(bin_data.iter().find(|b| b.0 == bin).is_some());
        }
        for (_bin, bin_start, bin_end, data) in bin_data.iter_mut() {
            let overlap_start = (*bin_start).max(interval_start);
            let overlap_end = (*bin_end).min(interval_end);

            for i in &mut data
                [((overlap_start - *bin_start) as usize)..((overlap_end - *bin_start) as usize)]
            {
                // If NAN, then 0.0 + 1.0, else i + 1.0
                *i = (*i).max(0.0) + 1.0;
            }
        }
    }
    while let Some(front) = bin_data.pop_front() {
        let bin = front.0;

        match summary {
            Summary::Min => {
                v[bin] = front
                    .3
                    .into_iter()
                    .reduce(|min, v| min.min(v))
                    .unwrap_or(missing);
            }
            Summary::Max => {
                v[bin] = front
                    .3
                    .into_iter()
                    .reduce(|max, v| max.max(v))
                    .unwrap_or(missing);
            }
            Summary::Mean => {
                v[bin] = front.3.into_iter().sum::<f64>() / bin_size;
            }
        }
    }
    Ok(())
}

fn to_entry_array_zoom<I: Iterator<Item = Result<ZoomRecord, BBIReadError>>>(
    start: i32,
    end: i32,
    iter: I,
    summary: Summary,
    bins: usize,
    missing: f64,
    mut v: ArrayViewMut<'_, f64, numpy::Ix1>,
) -> Result<(), BBIReadError> {
    assert_eq!(v.len(), bins);
    v.fill(missing);

    let mut bin_data: VecDeque<(usize, i32, i32, Vec<f64>)> = VecDeque::new();
    let bin_size = (end - start) as f64 / bins as f64;
    for interval in iter {
        let interval = interval?;
        let interval_start = (interval.start as i32).max(start) - start;
        let interval_end = (interval.end as i32).min(end) - start;
        let bin_start = ((interval_start as f64) / bin_size) as usize;
        let bin_end = ((interval_end as f64) / bin_size).ceil() as usize;

        while let Some(front) = bin_data.front_mut() {
            if front.0 < bin_start {
                let front = bin_data.pop_front().unwrap();
                let bin = front.0;

                match summary {
                    Summary::Min => {
                        v[bin] = front
                            .3
                            .into_iter()
                            .reduce(|min, v| min.min(v))
                            .unwrap_or(missing);
                    }
                    Summary::Max => {
                        v[bin] = front
                            .3
                            .into_iter()
                            .reduce(|max, v| max.max(v))
                            .unwrap_or(missing);
                    }
                    Summary::Mean => {
                        v[bin] = front.3.into_iter().sum::<f64>() / bin_size;
                    }
                }
            } else {
                break;
            }
        }
        while let Some(bin) = bin_data
            .back_mut()
            .map(|b| ((b.0 + 1) < bin_end).then(|| b.0 + 1))
            .unwrap_or(Some(bin_start))
        {
            let bin_start = ((bin as f64) * bin_size).ceil() as i32;
            let bin_end = (((bin + 1) as f64) * bin_size).ceil() as i32;

            bin_data.push_back((
                bin,
                bin_start,
                bin_end,
                vec![missing; (bin_end - bin_start) as usize],
            ));
        }
        for bin in bin_data.iter() {
            assert!(
                (bin_start..bin_end).contains(&bin.0),
                "{} not in {}..{}",
                bin.0,
                bin_start,
                bin_end
            );
        }
        for bin in bin_start..bin_end {
            assert!(bin_data.iter().find(|b| b.0 == bin).is_some());
        }
        for (_bin, bin_start, bin_end, data) in bin_data.iter_mut() {
            let overlap_start = (*bin_start).max(interval_start);
            let overlap_end = (*bin_end).min(interval_end);

            let mean = interval.summary.sum / (interval.end - interval.start) as f64;
            for i in &mut data
                [((overlap_start - *bin_start) as usize)..((overlap_end - *bin_start) as usize)]
            {
                *i = i.max(0.0);
                match summary {
                    Summary::Mean => *i += mean,
                    Summary::Min => *i = i.min(interval.summary.min_val),
                    Summary::Max => *i = i.max(interval.summary.max_val),
                }
            }
        }
    }
    while let Some(front) = bin_data.pop_front() {
        let bin = front.0;

        match summary {
            Summary::Min => {
                v[bin] = front
                    .3
                    .into_iter()
                    .reduce(|min, v| min.min(v))
                    .unwrap_or(missing);
            }
            Summary::Max => {
                v[bin] = front
                    .3
                    .into_iter()
                    .reduce(|max, v| max.max(v))
                    .unwrap_or(missing);
            }
            Summary::Mean => {
                v[bin] = front.3.into_iter().sum::<f64>() / bin_size;
            }
        }
    }
    Ok(())
}

/// This class is the interface for reading a bigWig or bigBed.
#[pyclass(module = "pybigtools")]
struct BBIRead {
    bbi: BBIReadRaw,
}

#[pymethods]
impl BBIRead {
    #[getter]
    fn is_bigwig(&self) -> bool {
        #[cfg(feature = "remote")]
        {
            matches!(
                self.bbi,
                BBIReadRaw::BigWigFile(_)
                    | BBIReadRaw::BigWigRemote(_)
                    | BBIReadRaw::BigWigFileLike(_)
            )
        }
        #[cfg(not(feature = "remote"))]
        {
            matches!(
                self.bbi,
                BBIReadRaw::BigWigFile(_) | BBIReadRaw::BigWigFileLike(_)
            )
        }
    }

    #[getter]
    fn is_bigbed(&self) -> bool {
        #[cfg(feature = "remote")]
        {
            matches!(
                self.bbi,
                BBIReadRaw::BigBedFile(_)
                    | BBIReadRaw::BigBedRemote(_)
                    | BBIReadRaw::BigBedFileLike(_)
            )
        }
        #[cfg(not(feature = "remote"))]
        {
            matches!(
                self.bbi,
                BBIReadRaw::BigBedFile(_) | BBIReadRaw::BigBedFileLike(_)
            )
        }
    }

    fn info(&mut self, py: Python<'_>) -> PyResult<PyObject> {
        let (info, summary) = match &mut self.bbi {
            BBIReadRaw::BigWigFile(b) => {
                let summary = b.get_summary()?;
                (b.info(), summary)
            }
            BBIReadRaw::BigWigRemote(b) => {
                let summary = b.get_summary()?;
                (b.info(), summary)
            }
            BBIReadRaw::BigWigFileLike(b) => {
                let summary = b.get_summary()?;
                (b.info(), summary)
            }
            BBIReadRaw::BigBedFile(b) => {
                let summary = b.get_summary()?;
                (b.info(), summary)
            }
            BBIReadRaw::BigBedRemote(b) => {
                let summary = b.get_summary()?;
                (b.info(), summary)
            }
            BBIReadRaw::BigBedFileLike(b) => {
                let summary = b.get_summary()?;
                (b.info(), summary)
            }
        };
        let var = (summary.sum_squares
            - (summary.sum * summary.sum) / summary.bases_covered as f64)
            / (summary.bases_covered as f64 - 1.0);
        let summary = [
            ("basesCovered", summary.bases_covered.to_object(py)),
            ("sum", summary.sum.to_object(py)),
            (
                "mean",
                (summary.sum as f64 / summary.bases_covered as f64).to_object(py),
            ),
            ("min", summary.min_val.to_object(py)),
            ("max", summary.max_val.to_object(py)),
            ("std", f64::sqrt(var).to_object(py)),
        ]
        .into_py_dict(py)
        .to_object(py);
        let info = [
            ("version", info.header.version.to_object(py)),
            ("isCompressed", info.header.is_compressed().to_object(py)),
            (
                "primaryDataSize",
                info.header.primary_data_size().to_object(py),
            ),
            ("zoomLevels", info.zoom_headers.len().to_object(py)),
            ("chromCount", info.chrom_info.len().to_object(py)),
            ("summary", summary),
        ]
        .into_py_dict(py)
        .to_object(py);
        Ok(info)
    }

    /// Returns the autosql of this bbi file.
    ///
    /// For bigBeds, this comes directly from the autosql stored in the file.
    /// For bigWigs, the autosql returned matches that of a bedGraph file.
    ///
    /// By default, the autosql is returned as a string. Passing `parse = true`
    /// returns instead a dictionary of the format:
    /// ```
    /// {
    ///   "name": <declared name>,
    ///   "comment": <declaration coment>,
    ///   "fields": [(<field name>, <field type>, <field comment>), ...],
    /// }
    /// ```
    fn sql(&mut self, py: Python, parse: Option<bool>) -> PyResult<PyObject> {
        let parse = parse.unwrap_or(false);
        pub const BEDGRAPH: &str = r#"
            table bedGraph
            "bedGraph file"
            (
                string chrom;        "Reference sequence chromosome or scaffold"
                uint   chromStart;   "Start position in chromosome"
                uint   chromEnd;     "End position in chromosome"
                float  value;        "Value for a given interval"
            )
            "#;
        let schema = match &mut self.bbi {
            BBIReadRaw::BigWigFile(_)
            | BBIReadRaw::BigWigRemote(_)
            | BBIReadRaw::BigWigFileLike(_) => BEDGRAPH.to_string(),
            BBIReadRaw::BigBedFile(b) => b.autosql().convert_err()?,
            BBIReadRaw::BigBedRemote(b) => b.autosql().convert_err()?,
            BBIReadRaw::BigBedFileLike(b) => b.autosql().convert_err()?,
        };
        let obj = if parse {
            let mut declarations = parse_autosql(&schema).map_err(|_| {
                PyErr::new::<exceptions::PyException, _>("Unable to parse autosql.")
            })?;
            if declarations.len() > 1 {
                return Err(PyErr::new::<exceptions::PyException, _>(
                    "Unexpected extra declarations.",
                ));
            }
            let declaration = declarations.pop();
            match declaration {
                None => PyDict::new(py).to_object(py),
                Some(d) => {
                    let fields = d
                        .fields
                        .iter()
                        .map(|f| (&f.name, f.field_type.to_string(), &f.comment))
                        .collect::<Vec<_>>()
                        .to_object(py);
                    [
                        ("name", d.name.name.to_object(py)),
                        ("comment", d.comment.to_object(py)),
                        ("fields", fields),
                    ]
                    .into_py_dict(py)
                    .to_object(py)
                }
            }
        } else {
            schema.to_object(py)
        };
        Ok(obj)
    }

    /// Returns the records of a given range on a chromosome.
    ///
    /// The result is an iterator of tuples. For bigWigs, these tuples are in
    /// the format (start: int, end: int, value: float). For bigBeds, these
    /// tuples are in the format (start: int, end: int, ...), where the "rest"
    /// fields are split by whitespace.
    ///
    /// Missing values in bigWigs will results in non-contiguous records.
    ///
    /// The chrom argument is the name of the chromosome.  
    /// The start and end arguments denote the range to get values for.
    ///  If end is not provided, it defaults to the length of the chromosome.
    ///  If start is not provided, it defaults to the beginning of the chromosome.
    fn records(
        &mut self,
        py: Python<'_>,
        chrom: String,
        start: Option<u32>,
        end: Option<u32>,
    ) -> PyResult<PyObject> {
        let (start, end) = start_end(&self.bbi, &chrom, start, end)?;
        match &self.bbi {
            BBIReadRaw::BigWigFile(b) => {
                let b = b.reopen()?;
                Ok(BigWigIntervalIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).unwrap()),
                }
                .into_py(py))
            }
            #[cfg(feature = "remote")]
            BBIReadRaw::BigWigRemote(b) => {
                let b = b.reopen()?;
                Ok(BigWigIntervalIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).unwrap()),
                }
                .into_py(py))
            }
            BBIReadRaw::BigWigFileLike(b) => {
                let b = b.reopen()?;
                Ok(BigWigIntervalIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).unwrap()),
                }
                .into_py(py))
            }
            BBIReadRaw::BigBedFile(b) => {
                let b = b.reopen()?;
                Ok(BigBedEntriesIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).unwrap()),
                }
                .into_py(py))
            }
            #[cfg(feature = "remote")]
            BBIReadRaw::BigBedRemote(b) => {
                let b = b.reopen()?;
                Ok(BigBedEntriesIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).unwrap()),
                }
                .into_py(py))
            }
            BBIReadRaw::BigBedFileLike(b) => {
                let b = b.reopen()?;
                Ok(BigBedEntriesIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).unwrap()),
                }
                .into_py(py))
            }
        }
    }

    /// Returns the values of a given range on a chromosome.
    ///
    /// For bigWigs, the result is an array of length (end - start).
    /// If a value does not exist in the bigwig for a specific base, it will be nan.
    ///
    /// For bigBeds, the returned array instead represents a pileup of the count
    /// of intervals overlapping each base.
    ///
    /// The chrom argument is the name of the chromosome.  
    /// The start and end arguments denote the range to get values for.  
    ///  If end is not provided, it defaults to the length of the chromosome.  
    ///  If start is not provided, it defaults to the beginning of the chromosome.  
    ///
    /// This returns a numpy array.
    #[pyo3(signature = (chrom, start, end, bins=None, summary="mean".to_string(), exact=false, missing=0.0, oob=f64::NAN, arr=None))]
    fn values(
        &mut self,
        py: Python<'_>,
        chrom: String,
        start: Option<i32>,
        end: Option<i32>,
        bins: Option<usize>,
        summary: String,
        exact: bool,
        missing: f64,
        oob: f64,
        arr: Option<PyObject>,
    ) -> PyResult<PyObject> {
        let summary = match summary.as_ref() {
            "mean" => Summary::Mean,
            "min" => Summary::Min,
            "max" => Summary::Max,
            _ => {
                return Err(PyErr::new::<exceptions::PyException, _>(format!(
                    "Unrecognized summary. Only `mean`, `min`, and `max` are allowed."
                )));
            }
        };
        match &mut self.bbi {
            BBIReadRaw::BigWigFile(b) => intervals_to_array(
                py, b, &chrom, start, end, bins, summary, exact, missing, oob, arr,
            ),
            #[cfg(feature = "remote")]
            BBIReadRaw::BigWigRemote(b) => intervals_to_array(
                py, b, &chrom, start, end, bins, summary, exact, missing, oob, arr,
            ),
            BBIReadRaw::BigWigFileLike(b) => intervals_to_array(
                py, b, &chrom, start, end, bins, summary, exact, missing, oob, arr,
            ),
            BBIReadRaw::BigBedFile(b) => entries_to_array(
                py, b, &chrom, start, end, bins, summary, exact, missing, oob, arr,
            ),
            #[cfg(feature = "remote")]
            BBIReadRaw::BigBedRemote(b) => entries_to_array(
                py, b, &chrom, start, end, bins, summary, exact, missing, oob, arr,
            ),
            BBIReadRaw::BigBedFileLike(b) => entries_to_array(
                py, b, &chrom, start, end, bins, summary, exact, missing, oob, arr,
            ),
        }
    }

    /// Returns the chromosomes in a bigwig, and their lengths.  
    ///
    /// The chroms argument can be either String or None.  
    ///  If it is None, then all chroms will be returned.  
    ///  If it is a String, then the length of that chromosome will be returned.  
    ///  If the chromosome doesn't exist, nothing will be returned.  
    fn chroms(&mut self, py: Python, chrom: Option<String>) -> Option<PyObject> {
        fn get_chrom_obj<B: bigtools::BBIRead>(
            b: &B,
            py: Python,
            chrom: Option<String>,
        ) -> Option<PyObject> {
            match chrom {
                Some(chrom) => b
                    .chroms()
                    .into_iter()
                    .find(|c| c.name == chrom)
                    .map(|c| c.length)
                    .map(|c| c.to_object(py)),
                None => Some(
                    b.chroms()
                        .into_iter()
                        .map(|c| (c.name.clone(), c.length))
                        .into_py_dict(py)
                        .into(),
                ),
            }
        }
        match &self.bbi {
            BBIReadRaw::BigWigFile(b) => get_chrom_obj(b, py, chrom),
            #[cfg(feature = "remote")]
            BBIReadRaw::BigWigRemote(b) => get_chrom_obj(b, py, chrom),
            BBIReadRaw::BigWigFileLike(b) => get_chrom_obj(b, py, chrom),
            BBIReadRaw::BigBedFile(b) => get_chrom_obj(b, py, chrom),
            #[cfg(feature = "remote")]
            BBIReadRaw::BigBedRemote(b) => get_chrom_obj(b, py, chrom),
            BBIReadRaw::BigBedFileLike(b) => get_chrom_obj(b, py, chrom),
        }
    }
}

/// This class is an iterator for `Values` from a bigWig.  
/// It returns only values that exist in the bigWig, skipping
/// any missing intervals.
#[pyclass(module = "pybigtools")]
struct BigWigIntervalIterator {
    iter: Box<dyn Iterator<Item = Result<Value, BBIReadError>> + Send>,
}

#[pymethods]
impl BigWigIntervalIterator {
    fn __iter__(slf: PyRefMut<Self>) -> PyResult<Py<BigWigIntervalIterator>> {
        Ok(slf.into())
    }

    fn __next__(mut slf: PyRefMut<Self>) -> PyResult<Option<ValueTuple>> {
        slf.iter
            .next()
            .transpose()
            .map(|o| o.map(|v| (v.start, v.end, v.value)))
            .convert_err()
    }
}

/// This class is an interator for the entries in a bigBed file.
#[pyclass(module = "pybigtools")]
struct BigBedEntriesIterator {
    iter: Box<dyn Iterator<Item = Result<BedEntry, BBIReadError>> + Send>,
}

#[pymethods]
impl BigBedEntriesIterator {
    fn __iter__(slf: PyRefMut<Self>) -> PyResult<Py<BigBedEntriesIterator>> {
        Ok(slf.into())
    }

    fn __next__(mut slf: PyRefMut<Self>) -> PyResult<Option<PyObject>> {
        let py = slf.py();
        let next = match slf.iter.next() {
            Some(n) => n.convert_err()?,
            None => return Ok(None),
        };
        let elements: Vec<_> = [next.start.to_object(py), next.end.to_object(py)]
            .into_iter()
            .chain(next.rest.split_whitespace().map(|o| o.to_object(py)))
            .collect();
        Ok(Some(
            PyTuple::new::<PyObject, _>(py, elements.into_iter()).to_object(py),
        ))
    }
}

/// This class is the interface for writing a bigWig.
#[pyclass(module = "pybigtools")]
struct BigWigWrite {
    bigwig: Option<BigWigWriteRaw>,
}

#[pymethods]
impl BigWigWrite {
    /// Writes the values passsed to the bigwig file.
    /// The underlying file will be closed automatically when the function completes (and no other operations will be able to be performed).  
    ///
    /// The chroms argument should be a dictionary with keys as chromosome names and values as their length.  
    /// The vals argument should be an iterable with values (String, int, int, float) that represents each value to write in the format (chromosome, start, end, value).  
    fn write(&mut self, py: Python, chroms: &PyDict, vals: Py<PyAny>) -> PyResult<()> {
        let runtime = runtime::Builder::new_multi_thread()
            .worker_threads(
                std::thread::available_parallelism()
                    .map(|c| c.into())
                    .unwrap_or(1),
            )
            .build()
            .expect("Unable to create thread pool.");

        let bigwig = self
            .bigwig
            .take()
            .ok_or_else(|| PyErr::new::<exceptions::PyException, _>("File already closed."))?;
        let chrom_map = chroms
            .into_iter()
            .map(|(key, val)| {
                let chrom: String = key.downcast::<PyString>()?.to_str().unwrap().to_owned();
                let length: u32 = val.downcast::<PyInt>()?.to_object(py).extract(py).unwrap();
                Ok((chrom, length))
            })
            .collect::<Result<std::collections::HashMap<String, u32>, pyo3::PyDowncastError>>()?;
        struct IterError(String);
        struct Iter {
            inner: PyObject,
        }
        impl Iterator for Iter {
            type Item = Result<(String, Value), IterError>;
            fn next(&mut self) -> Option<Self::Item> {
                // We have to reacquire the gil for each iteration
                Python::with_gil(|py| {
                    let mut iter: &PyIterator = match self.inner.downcast(py) {
                        Ok(o) => o,
                        Err(_) => {
                            return Some(Err(IterError(format!(
                                "Passed value for `val` is not iterable."
                            ))))
                        }
                    };
                    let next: Result<(String, Value), pyo3::PyDowncastError> = match iter.next()? {
                        Err(e) => {
                            e.print(py);
                            return Some(Err(IterError(format!(
                                "An error occurred while iterating."
                            ))));
                        }
                        Ok(n) => {
                            // TODO: try block or separate function
                            (|| {
                                let tuple = n.downcast::<PyTuple>()?;
                                assert!(tuple.len() == 4);
                                let chrom: String = tuple
                                    .get_item(0)
                                    .unwrap()
                                    .downcast::<PyString>()?
                                    .to_str()
                                    .unwrap()
                                    .to_owned();
                                let start: u32 = tuple
                                    .get_item(1)
                                    .unwrap()
                                    .downcast::<PyInt>()?
                                    .to_object(py)
                                    .extract(py)
                                    .unwrap();
                                let end: u32 = tuple
                                    .get_item(2)
                                    .unwrap()
                                    .downcast::<PyInt>()?
                                    .to_object(py)
                                    .extract(py)
                                    .unwrap();
                                let value: f32 = tuple
                                    .get_item(3)
                                    .unwrap()
                                    .downcast::<PyFloat>()?
                                    .to_object(py)
                                    .extract(py)
                                    .unwrap();
                                Ok((chrom, Value { start, end, value }))
                            })()
                        }
                    };
                    let ret = match next {
                        Err(_) => Err(IterError("Invalid iterator value. Must a tuple of type (String, int, int, float)".to_string())),
                        Ok(n) => Ok(n),
                    };
                    Some(ret)
                })
            }
        }
        py.allow_threads(|| {
            let iter = Python::with_gil(|py| {
                let inner_obj: PyObject = vals.into_py(py);
                match PyIterator::from_object(py, &inner_obj) {
                    Ok(iter) => Ok(iter.to_object(py)),
                    Err(_) => Err(PyTypeError::new_err(
                        "Passed value for `val` is not iterable.",
                    )),
                }
            })?;
            let vals_iter_raw = Iter { inner: iter }.map(|v| match v {
                Err(e) => Err(io::Error::new(io::ErrorKind::Other, format!("{}", e.0))),
                Ok(v) => Ok(v),
            });
            let vals_iter = BedParser::wrap_iter(vals_iter_raw);
            let chsi = BedParserStreamingIterator::new(vals_iter, true);
            match bigwig.write(chrom_map, chsi, runtime) {
                Err(e) => println!("{}", e),
                Ok(_) => {}
            }
            Ok(())
        })
    }

    /// close()
    /// --
    ///
    /// Manually closed the file.
    /// No other operations will be allowed after it is closed. This is done automatically after write is performed.
    fn close(&mut self) -> PyResult<()> {
        self.bigwig.take();
        Ok(())
    }
}

/// This class is the interface for writing to a bigBed.
#[pyclass(module = "pybigtools")]
struct BigBedWrite {
    bigbed: Option<BigBedWriteRaw>,
}

#[pymethods]
impl BigBedWrite {
    /// Writes the values passsed to the bigwig file. The underlying file will be closed automatically when the function completes (and no other operations will be able to be performed).  
    /// The chroms argument should be a dictionary with keys as chromosome names and values as their length.  
    /// The vals argument should be an iterable with values (String, int, int, String) that represents each value to write in the format (chromosome, start, end, rest)  
    ///   The rest String should consist of tab-delimited fields.  
    fn write(&mut self, py: Python, chroms: &PyDict, vals: Py<PyAny>) -> PyResult<()> {
        let runtime = runtime::Builder::new_multi_thread()
            .worker_threads(
                std::thread::available_parallelism()
                    .map(|c| c.into())
                    .unwrap_or(1),
            )
            .build()
            .expect("Unable to create thread pool.");

        let bigbed = self
            .bigbed
            .take()
            .ok_or_else(|| PyErr::new::<exceptions::PyException, _>("File already closed."))?;
        let chrom_map = chroms
            .into_iter()
            .map(|(key, val)| {
                let chrom: String = key.downcast::<PyString>()?.to_str().unwrap().to_owned();
                let length: u32 = val.downcast::<PyInt>()?.to_object(py).extract(py).unwrap();
                Ok((chrom, length))
            })
            .collect::<Result<std::collections::HashMap<String, u32>, pyo3::PyDowncastError>>()?;
        struct IterError(String);
        struct Iter {
            inner: PyObject,
        }
        impl Iterator for Iter {
            type Item = Result<(String, BedEntry), IterError>;
            fn next(&mut self) -> Option<Self::Item> {
                // We have to reacquire the gil for each iteration
                Python::with_gil(|py| {
                    let mut iter: &PyIterator = match self.inner.downcast(py) {
                        Ok(o) => o,
                        Err(_) => {
                            return Some(Err(IterError(format!(
                                "Passed value for `val` is not iterable."
                            ))))
                        }
                    };
                    let next: Result<(String, BedEntry), pyo3::PyDowncastError> =
                        match iter.next()? {
                            Err(e) => {
                                e.print(py);
                                return Some(Err(IterError(format!(
                                    "An error occurred while iterating."
                                ))));
                            }
                            Ok(n) => {
                                // TODO: try block or separate function
                                (|| {
                                    let tuple = n.downcast::<PyTuple>()?;
                                    assert!(tuple.len() == 4);
                                    let chrom: String = tuple
                                        .get_item(0)
                                        .unwrap()
                                        .downcast::<PyString>()?
                                        .to_str()
                                        .unwrap()
                                        .to_owned();
                                    let start: u32 = tuple
                                        .get_item(1)
                                        .unwrap()
                                        .downcast::<PyInt>()?
                                        .to_object(py)
                                        .extract(py)
                                        .unwrap();
                                    let end: u32 = tuple
                                        .get_item(2)
                                        .unwrap()
                                        .downcast::<PyInt>()?
                                        .to_object(py)
                                        .extract(py)
                                        .unwrap();
                                    let rest: String = tuple
                                        .get_item(0)
                                        .unwrap()
                                        .downcast::<PyString>()?
                                        .to_str()
                                        .unwrap()
                                        .to_owned();
                                    Ok((chrom, BedEntry { start, end, rest }))
                                })()
                            }
                        };
                    let ret = match next {
                        Err(_) => Err(IterError("Invalid iterator value. Must a tuple of type (String, int, int, String)".to_string())),
                        Ok(n) => Ok(n),
                    };
                    Some(ret)
                })
            }
        }
        py.allow_threads(|| {
            let iter = Python::with_gil(|py| {
                let inner_obj: PyObject = vals.into_py(py);
                match PyIterator::from_object(py, &inner_obj) {
                    Ok(iter) => Ok(iter.to_object(py)),
                    Err(_) => Err(PyTypeError::new_err(
                        "Passed value for `val` is not iterable.",
                    )),
                }
            })?;
            let vals_iter_raw = Iter { inner: iter }.map(|v| match v {
                Err(e) => Err(io::Error::new(io::ErrorKind::Other, format!("{}", e.0))),
                Ok(v) => Ok(v),
            });
            let vals_iter = BedParser::wrap_iter(vals_iter_raw);
            let chsi = BedParserStreamingIterator::new(vals_iter, true);
            match bigbed.write(chrom_map, chsi, runtime) {
                Err(e) => {
                    println!("{}", e)
                }
                Ok(_) => {}
            }
            Ok(())
        })
    }

    /// Manually closed the file. No other operations will be allowed after it is closed. This is done automatically after write is performed.
    fn close(&mut self) -> PyResult<()> {
        self.bigbed.take();
        Ok(())
    }
}

enum BigWigAverageOverBedEntriesIteratorRet {
    Single(f64),
    WithName((String, f64)),
}

impl IntoPy<PyObject> for BigWigAverageOverBedEntriesIteratorRet {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            BigWigAverageOverBedEntriesIteratorRet::Single(v) => v.into_py(py),
            BigWigAverageOverBedEntriesIteratorRet::WithName(v) => v.into_py(py),
        }
    }
}

/// This class is an interator for the entries of bigWigAverageOverBed
#[pyclass(module = "pybigtools")]
struct BigWigAverageOverBedEntriesIterator {
    iter: Box<
        dyn Iterator<Item = Result<BigWigAverageOverBedEntry, BigWigAverageOverBedError>> + Send,
    >,
    usename: bool,
}

#[pymethods]
impl BigWigAverageOverBedEntriesIterator {
    fn __iter__(slf: PyRefMut<Self>) -> PyResult<Py<BigWigAverageOverBedEntriesIterator>> {
        Ok(slf.into())
    }

    fn __next__(
        mut slf: PyRefMut<Self>,
    ) -> PyResult<Option<BigWigAverageOverBedEntriesIteratorRet>> {
        slf.iter
            .next()
            .transpose()
            .map(|o| {
                o.map(|v| {
                    if slf.usename {
                        BigWigAverageOverBedEntriesIteratorRet::WithName((v.name, v.mean))
                    } else {
                        BigWigAverageOverBedEntriesIteratorRet::Single(v.mean)
                    }
                })
            })
            .map_err(|e| PyErr::new::<exceptions::PyException, _>(format!("{}", e)))
    }
}

/// This is the entrypoint for working with bigWigs or bigBeds.
///
/// The first argument can take one of three values:
/// 1) A path to a file as a string
/// 2) An http url for a remote file as a string
/// 3) A file-like object with a `read` and `seek` method
///
/// When writing, only a file path can be used.
///
/// The mode is either "w", "r", or None. "r" or None will open a
/// bigWig or bigBed for reading (but will not allow writing).
/// "w" Will open a bigWig/bigBed for writing (but will not allow reading).
///
/// If passing a file-like object, simultaneous reading of different intervals
/// is not supported and may result in incorrect behavior.
///
/// This returns one of the following:  
/// - `BigWigWrite`
/// - `BigBedWrite`
/// - `BBIRead`
#[pyfunction]
fn open(py: Python, path_url_or_file_like: PyObject, mode: Option<String>) -> PyResult<PyObject> {
    let iswrite = match &mode {
        Some(mode) if mode == "w" => true,
        Some(mode) if mode == "r" => false,
        None => false,
        Some(mode) => {
            return Err(PyErr::new::<exceptions::PyException, _>(format!(
                "Invalid mode: `{}`",
                mode
            )));
        }
    };

    // If string, might be path or url like
    if let Ok(string_ref) = path_url_or_file_like.downcast::<PyString>(py) {
        return open_path_or_url(py, string_ref.to_str().unwrap().to_owned(), iswrite);
    }

    if iswrite {
        return Err(PyErr::new::<exceptions::PyException, _>(format!(
            "Writing only supports file names",
        )));
    }
    let file_like = match PyFileLikeObject::new(path_url_or_file_like, true, false, true) {
        Ok(file_like) => file_like,
        Err(_) => return Err(PyErr::new::<exceptions::PyException, _>(format!(
            "Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.",
        ))),
    };
    let read = match BigWigReadRaw::open(file_like.clone()) {
        Ok(bwr) => BBIRead {
            bbi: BBIReadRaw::BigWigFileLike(bwr.cached()),
        }
        .into_py(py),
        Err(_) => match BigBedReadRaw::open(file_like) {
            Ok(bbr) => BBIRead {
                bbi: BBIReadRaw::BigBedFileLike(bbr.cached()),
            }
            .into_py(py),
            Err(e) => {
                return Err(PyErr::new::<exceptions::PyException, _>(format!(
                    "File-like object is not a bigWig or bigBed. Or there was just a problem reading: {e}",
                )))
            }
        },
    };
    Ok(read)
}

fn open_path_or_url(
    py: Python,
    path_url_or_file_like: String,
    iswrite: bool,
) -> PyResult<PyObject> {
    let extension = match &Path::new(&path_url_or_file_like)
        .extension()
        .map(|e| e.to_string_lossy())
    {
        Some(e) => e.to_string(),
        None => {
            return Err(PyErr::new::<exceptions::PyException, _>(format!("Invalid file type. Must be either a bigWig (.bigWig, .bw) or bigBed (.bigBed, .bb).")));
        }
    };
    let res = if iswrite {
        match Url::parse(&path_url_or_file_like) {
            Ok(_) => {
                return Err(PyErr::new::<exceptions::PyException, _>(format!(
                    "Invalid file path. Writing does not support urls."
                )))
            }
            Err(_) => {}
        }
        match extension.as_ref() {
            "bw" | "bigWig" | "bigwig" => BigWigWrite {
                bigwig: Some(BigWigWriteRaw::create_file(path_url_or_file_like)),
            }
            .into_py(py),
            "bb" | "bigBed" | "bigbed" => BigBedWrite {
                bigbed: Some(BigBedWriteRaw::create_file(path_url_or_file_like)),
            }
            .into_py(py),
            _ => {
                return Err(PyErr::new::<exceptions::PyException, _>(format!("Invalid file type. Must be either a bigWig (.bigWig, .bw) or bigBed (.bigBed, .bb).")));
            }
        }
    } else {
        let isfile = Path::new(&path_url_or_file_like).exists();
        if !isfile {
            match Url::parse(&path_url_or_file_like) {
                Ok(_) => {}
                Err(_) => {
                    return Err(PyErr::new::<exceptions::PyException, _>(format!(
                        "Invalid file path. The file does not exists and it is not a url."
                    )))
                }
            }
        }
        match extension.as_ref() {
            "bw" | "bigWig" | "bigwig" => {
                if isfile {
                    match BigWigReadRaw::open_file(&path_url_or_file_like) {
                        Ok(bwr) => BBIRead {
                            bbi: BBIReadRaw::BigWigFile(bwr.cached()),
                        }
                        .into_py(py),
                        Err(_) => {
                            return Err(PyErr::new::<exceptions::PyException, _>(format!(
                                "Error opening bigWig."
                            )))
                        }
                    }
                } else {
                    #[cfg(feature = "remote")]
                    match BigWigReadRaw::open(RemoteFile::new(&path_url_or_file_like)) {
                        Ok(bwr) => BBIRead {
                            bbi: BBIReadRaw::BigWigRemote(bwr.cached()),
                        }
                        .into_py(py),
                        Err(_) => {
                            return Err(PyErr::new::<exceptions::PyException, _>(format!(
                                "Error opening bigWig."
                            )))
                        }
                    }

                    #[cfg(not(feature = "remote"))]
                    {
                        return Err(PyErr::new::<exceptions::PyException, _>(format!(
                            "Builtin support for remote files is not supported on this platform."
                        )));
                    }
                }
            }
            "bb" | "bigBed" | "bigbed" => {
                if isfile {
                    match BigBedReadRaw::open_file(&path_url_or_file_like) {
                        Ok(bwr) => BBIRead {
                            bbi: BBIReadRaw::BigBedFile(bwr.cached()),
                        }
                        .into_py(py),
                        Err(_) => {
                            return Err(PyErr::new::<exceptions::PyException, _>(format!(
                                "Error opening bigBed."
                            )))
                        }
                    }
                } else {
                    #[cfg(feature = "remote")]
                    match BigBedReadRaw::open(RemoteFile::new(&path_url_or_file_like)) {
                        Ok(bwr) => BBIRead {
                            bbi: BBIReadRaw::BigBedRemote(bwr.cached()),
                        }
                        .into_py(py),
                        Err(_) => {
                            return Err(PyErr::new::<exceptions::PyException, _>(format!(
                                "Error opening bigBed."
                            )))
                        }
                    }

                    #[cfg(not(feature = "remote"))]
                    {
                        return Err(PyErr::new::<exceptions::PyException, _>(format!(
                            "Builtin support for remote files is not supported on this platform."
                        )));
                    }
                }
            }
            _ => {
                return Err(PyErr::new::<exceptions::PyException, _>(format!("Invalid file type. Must be either a bigWig (.bigWig, .bw) or bigBed (.bigBed, .bb).")));
            }
        }
    };
    Ok(res)
}

/// Gets the average values from a bigWig over the entries of a bed file.
///
/// Parameters:
///     bigWig (str): The path to the bigWig.
///     bed (str): The path to the bed.
///     names (None, bool, or int):  
///         If `None`, then each return value will be a single `float`,
///             the average value over an interval in the bed file.  
///         If `True`, then each return value will be a tuple of the value of column 4
///             and the average value over the interval with that name in the bed file.  
///         If `False`, then each return value will be a tuple of the interval in the format
///             `{chrom}:{start}-{end}` and the average value over that interval.  
///         If `0`, then each return value will match as if `False` was passed.  
///         If a `1+`, then each return value will be a tuple of the value of column of this
///             parameter (1-based) and the average value over the interval.  
///
/// Returns:
///     This returns a generator of values. (Therefore, to save to a list, do `list(bigWigAverageOverBed(...))`)  
///     If no name is specified (see the `names` parameter above), then returns a generator of `float`s.  
///     If a name column is specified (see above), then returns a generator of tuples `({name}, {average})`  
#[pyfunction]
fn bigWigAverageOverBed(
    py: Python,
    bigwig: String,
    bed: String,
    names: Option<PyObject>,
) -> PyResult<PyObject> {
    let extension = match &Path::new(&bigwig).extension().map(|e| e.to_string_lossy()) {
        Some(e) => e.to_string(),
        None => {
            return Err(PyErr::new::<exceptions::PyException, _>(format!(
                "Invalid file type. Must be a bigWig (.bigWig, .bw)."
            )));
        }
    };
    let isfile = std::path::Path::new(&bigwig).exists();
    if !isfile {
        match Url::parse(&bigwig) {
            Ok(_) => {}
            Err(_) => {
                return Err(PyErr::new::<exceptions::PyException, _>(format!(
                    "Invalid file path. The file does not exists and it is not a url."
                )))
            }
        }
    }
    let (name, usename) = {
        match names {
            Some(names) => match names.extract::<bool>(py) {
                Ok(b) => {
                    if b {
                        (Name::Column(3), true)
                    } else {
                        (Name::None, true)
                    }
                }
                Err(_) => match names.extract::<isize>(py) {
                    Ok(col) => match col {
                        0 => (Name::None, true),
                        1.. => (Name::Column((col - 1) as usize), true),
                        _ => {
                            return Err(PyErr::new::<exceptions::PyException, _>(
                                "Invalid names argument. Must be >= 0.",
                            ));
                        }
                    },
                    Err(_) => {
                        return Err(PyErr::new::<exceptions::PyException, _>("Invalid names argument. Should be either `None`, a `bool`, or an `int`"));
                    }
                },
            },
            None => (Name::None, false),
        }
    };
    let res = match extension.as_ref() {
        "bw" | "bigWig" | "bigwig" => {
            if isfile {
                let read = BigWigReadRaw::open_file(&bigwig)
                    .map_err(|_| {
                        PyErr::new::<exceptions::PyException, _>(format!("Error opening bigWig."))
                    })?
                    .cached();
                let bedin = BufReader::new(File::open(bed)?);
                let iter = Box::new(bigwig_average_over_bed(bedin, read, name));

                BigWigAverageOverBedEntriesIterator { iter, usename }.into_py(py)
            } else {
                #[cfg(feature = "remote")]
                {
                    let read = BigWigReadRaw::open(RemoteFile::new(&bigwig))
                        .map_err(|_| {
                            PyErr::new::<exceptions::PyException, _>(format!(
                                "Error opening bigBed."
                            ))
                        })?
                        .cached();
                    let bedin = BufReader::new(File::open(bed)?);
                    let iter = Box::new(bigwig_average_over_bed(bedin, read, name));

                    BigWigAverageOverBedEntriesIterator { iter, usename }.into_py(py)
                }

                #[cfg(not(feature = "remote"))]
                {
                    return Err(PyErr::new::<exceptions::PyException, _>(format!(
                        "Builtin support for remote files is not available on this platform."
                    )));
                }
            }
        }
        _ => {
            return Err(PyErr::new::<exceptions::PyException, _>(format!(
                "Invalid file type. Must be a bigWig (.bigWig, .bw)."
            )));
        }
    };

    Ok(res)
}

/// The base module for opening a bigWig or bigBed. The only defined function is `open`.
#[pymodule]
fn pybigtools(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_wrapped(wrap_pyfunction!(open))?;
    m.add_wrapped(wrap_pyfunction!(bigWigAverageOverBed))?;

    m.add_class::<BigWigWrite>()?;
    m.add_class::<BigBedWrite>()?;
    m.add_class::<BBIRead>()?;

    m.add_class::<BigWigIntervalIterator>()?;
    m.add_class::<BigBedEntriesIterator>()?;

    Ok(())
}

#[cfg(test)]
mod test {
    use bigtools::BedEntry;
    use numpy::ndarray::Array;

    use crate::to_entry_array_bins;

    #[test]
    fn test_to_entry_array_bins() {
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
        to_entry_array_bins(
            10,
            20,
            entries.into_iter().map(|v| Ok(v)),
            crate::Summary::Mean,
            1,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [f64::NAN].into_iter().collect::<Vec<_>>());

        let entries = [];
        let mut arr = Array::from(vec![f64::NAN]);
        to_entry_array_bins(
            10,
            20,
            entries.into_iter().map(|v| Ok(v)),
            crate::Summary::Mean,
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
        to_entry_array_bins(
            10,
            20,
            entries.into_iter().map(|v| Ok(v)),
            crate::Summary::Mean,
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
        to_entry_array_bins(
            10,
            20,
            entries.into_iter().map(|v| Ok(v)),
            crate::Summary::Mean,
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
        to_entry_array_bins(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            crate::Summary::Mean,
            2,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [1.0, 2.0].into_iter().collect::<Vec<_>>());
        let mut arr = Array::from(vec![f64::NAN]);
        to_entry_array_bins(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            crate::Summary::Min,
            1,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [1.0].into_iter().collect::<Vec<_>>());
        let mut arr = Array::from(vec![f64::NAN]);
        to_entry_array_bins(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            crate::Summary::Max,
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
        to_entry_array_bins(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            crate::Summary::Mean,
            2,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [1.0, 3.0].into_iter().collect::<Vec<_>>());
        let mut arr = Array::from(vec![f64::NAN]);
        to_entry_array_bins(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            crate::Summary::Min,
            1,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [1.0].into_iter().collect::<Vec<_>>());
        let mut arr = Array::from(vec![f64::NAN]);
        to_entry_array_bins(
            10,
            20,
            entries.clone().into_iter().map(|v| Ok(v)),
            crate::Summary::Max,
            1,
            f64::NAN,
            arr.view_mut(),
        )
        .unwrap();
        let res = arr.to_vec();
        assert_equal(res, [3.0].into_iter().collect::<Vec<_>>());
    }
}
