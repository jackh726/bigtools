#![allow(non_snake_case)]

use std::fs::File;
use std::io::BufReader;
use std::ops::Deref;

use bigtools::bed::autosql::parse::parse_autosql;
use bigtools::utils::misc::{
    bigwig_average_over_bed, BigWigAverageOverBedEntry, BigWigAverageOverBedError, Name,
};
use bigtools::utils::reopen::Reopen;
use bigtools::{BBIReadError as _BBIReadError, BedEntry, Value, ZoomRecord};
use pyo3::exceptions::{self, PyKeyError};
use pyo3::types::{IntoPyDict, PyAny, PyDict, PyList, PyString, PyTuple};
use pyo3::{prelude::*, PyTraverseError, PyVisit};

use crate::errors::{BBIFileClosed, BBIReadError, ConvertResult};
use crate::utils::{entries_to_array, intervals_to_array, start_end, BBIReadRaw, Summary};

type ValueTuple = (u32, u32, f32);

enum BigWigAverageOverBedStatistics {
    Size,
    Bases,
    Sum,
    Mean0,
    Mean,
    Min,
    Max,
}

impl BigWigAverageOverBedStatistics {
    fn from_str(val: &str) -> Option<BigWigAverageOverBedStatistics> {
        Some(match val {
            "size" => BigWigAverageOverBedStatistics::Size,
            "bases" => BigWigAverageOverBedStatistics::Bases,
            "sum" => BigWigAverageOverBedStatistics::Sum,
            "mean0" => BigWigAverageOverBedStatistics::Mean0,
            "mean" => BigWigAverageOverBedStatistics::Mean,
            "min" => BigWigAverageOverBedStatistics::Min,
            "max" => BigWigAverageOverBedStatistics::Max,
            _ => {
                return None;
            }
        })
    }
}

/// This class is an interator for the entries of bigWigAverageOverBed
#[pyclass(module = "pybigtools")]
struct BigWigAverageOverBedEntriesIterator {
    iter: Box<
        dyn Iterator<Item = Result<(String, BigWigAverageOverBedEntry), BigWigAverageOverBedError>>
            + Send,
    >,
    usename: bool,
    stats: Option<Vec<BigWigAverageOverBedStatistics>>,
    summary_statistics: PyObject,
}

#[pymethods]
impl BigWigAverageOverBedEntriesIterator {
    fn __iter__(slf: PyRefMut<Self>) -> PyResult<Py<BigWigAverageOverBedEntriesIterator>> {
        Ok(slf.into())
    }

    fn __next__(mut slf: PyRefMut<Self>) -> PyResult<Option<PyObject>> {
        let v = slf
            .iter
            .next()
            .transpose()
            .map_err(|e| PyErr::new::<exceptions::PyException, _>(format!("{}", e)))?;

        let Some((name, v)) = v else {
            return Ok(None);
        };
        let stats = match &slf.stats {
            Some(stats) => {
                let mut ret = Vec::with_capacity(stats.len());
                for stat in stats.iter() {
                    match stat {
                        BigWigAverageOverBedStatistics::Size => {
                            ret.push(v.size.to_object(slf.py()))
                        }
                        BigWigAverageOverBedStatistics::Bases => {
                            ret.push(v.bases.to_object(slf.py()))
                        }
                        BigWigAverageOverBedStatistics::Sum => ret.push(v.sum.to_object(slf.py())),
                        BigWigAverageOverBedStatistics::Mean0 => {
                            ret.push(v.mean0.to_object(slf.py()))
                        }
                        BigWigAverageOverBedStatistics::Mean => {
                            ret.push(v.mean.to_object(slf.py()))
                        }
                        BigWigAverageOverBedStatistics::Min => ret.push(v.min.to_object(slf.py())),
                        BigWigAverageOverBedStatistics::Max => ret.push(v.max.to_object(slf.py())),
                    }
                }
                if ret.len() == 1 {
                    ret[0].to_object(slf.py())
                } else {
                    PyTuple::new_bound(slf.py(), ret).to_object(slf.py())
                }
            }
            None => {
                let val = slf.summary_statistics.call_bound(
                    slf.py(),
                    (v.size, v.bases, v.sum, v.mean0, v.mean, v.min, v.max),
                    None,
                )?;
                val.to_object(slf.py())
            }
        };

        match slf.usename {
            true => Ok(Some((name, stats).to_object(slf.py()))),
            false => Ok(Some(stats)),
        }
    }
}

/// Interface for reading a BigWig or BigBed file.
#[pyclass(module = "pybigtools")]
pub struct BBIRead {
    pub bbi: BBIReadRaw,
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

    /// Return a dict of information about the BBI file.
    fn info(&mut self, py: Python<'_>) -> PyResult<PyObject> {
        let (info, summary) = match &mut self.bbi {
            BBIReadRaw::Closed => return Err(BBIFileClosed::new_err("File is closed.")),
            BBIReadRaw::BigWigFile(b) => {
                let summary = b.get_summary()?;
                (b.info(), summary)
            }
            #[cfg(feature = "remote")]
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
            #[cfg(feature = "remote")]
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
        .into_py_dict_bound(py)
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
        .into_py_dict_bound(py)
        .to_object(py);
        Ok(info)
    }

    /// Return a list of sizes in bases of the summary intervals used in each
    /// of the zoom levels (i.e. reduction levels) of the BBI file.
    fn zooms(&self) -> PyResult<Vec<u32>> {
        let zooms = match &self.bbi {
            BBIReadRaw::Closed => return Err(BBIFileClosed::new_err("File is closed.")),
            BBIReadRaw::BigWigFile(b) => &b.info().zoom_headers,
            #[cfg(feature = "remote")]
            BBIReadRaw::BigWigRemote(b) => &b.info().zoom_headers,
            BBIReadRaw::BigWigFileLike(b) => &b.info().zoom_headers,
            BBIReadRaw::BigBedFile(b) => &b.info().zoom_headers,
            #[cfg(feature = "remote")]
            BBIReadRaw::BigBedRemote(b) => &b.info().zoom_headers,
            BBIReadRaw::BigBedFileLike(b) => &b.info().zoom_headers,
        };
        Ok(zooms.iter().map(|z| z.reduction_level).collect())
    }

    /// Return the autoSql schema definition of this BBI file.
    ///
    /// For BigBeds, this schema comes directly from the autoSql string stored
    /// in the file. For BigWigs, the schema generated describes a bedGraph
    /// file.
    ///
    /// Parameters
    /// ----------
    /// parse : bool, optional [default: False]
    ///     If True, return the schema as a dictionary. If False, return the
    ///     schema as a string. Default is False.
    ///
    /// Returns
    /// -------
    /// schema : str or dict
    ///     The autoSql schema of the BBI file. If `parse` is True, the schema
    ///     is returned as a dictionary of the format:
    ///
    ///     ```
    ///     {
    ///         "name": <declared name>,
    ///         "comment": <declaration coment>,
    ///         "fields": [(<field name>, <field type>, <field comment>), ...],
    ///     }
    ///     ```
    ///
    /// See Also
    /// --------
    /// is_bigwig : Check if the BBI file is a bigWig.
    /// is_bigbed : Check if the BBI file is a bigBed.
    /// info : Get information about the BBI file.
    /// zooms : Get the zoom levels of the BBI file.
    #[pyo3(signature = (parse = false))]
    fn sql(&mut self, py: Python, parse: bool) -> PyResult<PyObject> {
        pub const BEDGRAPH: &str = r#"table bedGraph
"bedGraph file"
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
    float  value;        "Value for a given interval"
)"#;
        let schema = match &mut self.bbi {
            BBIReadRaw::Closed => return Err(BBIFileClosed::new_err("File is closed.")),
            BBIReadRaw::BigWigFile(_) | BBIReadRaw::BigWigFileLike(_) => BEDGRAPH.to_string(),
            #[cfg(feature = "remote")]
            BBIReadRaw::BigWigRemote(_) => BEDGRAPH.to_string(),
            BBIReadRaw::BigBedFile(b) => b.autosql().convert_err()?.unwrap_or(String::new()),
            #[cfg(feature = "remote")]
            BBIReadRaw::BigBedRemote(b) => b.autosql().convert_err()?.unwrap_or(String::new()),
            BBIReadRaw::BigBedFileLike(b) => b.autosql().convert_err()?.unwrap_or(String::new()),
        };
        let obj = if parse {
            let mut declarations = parse_autosql(&schema)
                .map_err(|_| PyErr::new::<BBIReadError, _>("Unable to parse autosql."))?;
            if declarations.len() > 1 {
                return Err(PyErr::new::<BBIReadError, _>(
                    "Unexpected extra declarations.",
                ));
            }
            let declaration = declarations.pop();
            match declaration {
                None => PyDict::new_bound(py).to_object(py),
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
                    .into_py_dict_bound(py)
                    .to_object(py)
                }
            }
        } else {
            schema.to_object(py)
        };
        Ok(obj)
    }

    /// Return the records of a given range on a chromosome.
    ///
    /// The result is an iterator of tuples. For BigWigs, these tuples are in
    /// the format (start: int, end: int, value: float). For BigBeds, these
    /// tuples are in the format (start: int, end: int, ...), where the "rest"
    /// fields are split by whitespace.
    ///
    /// Parameters
    /// ----------
    /// chrom : str
    ///     Name of the chromosome.
    /// start, end : int, optional
    ///     The range to get values for. If end is not provided, it defaults to
    ///     the length of the chromosome. If start is not provided, it defaults
    ///     to the beginning of the chromosome.
    ///
    /// Returns
    /// -------
    /// Iterator[tuple[int, int, float] or tuple[int, int, ...]]
    ///     An iterator of tuples in the format (start: int, end: int, value:
    ///     float) for BigWigs, or (start: int, end: int, *rest) for BigBeds.
    ///
    /// Notes
    /// -----
    /// Missing values in BigWigs will results in non-contiguous records.
    ///
    /// See Also
    /// --------
    /// zoom_records : Get the zoom records of a given range on a chromosome.
    /// values : Get the values of a given range on a chromosome.
    #[pyo3(signature = (chrom, start=None, end=None))]
    fn records(
        &mut self,
        py: Python<'_>,
        chrom: String,
        start: Option<i32>,
        end: Option<i32>,
    ) -> PyResult<PyObject> {
        let (start, end) = start_end(&self.bbi, &chrom, start, end)?;
        match &self.bbi {
            BBIReadRaw::Closed => return Err(BBIFileClosed::new_err("File is closed.")),
            BBIReadRaw::BigWigFile(b) => {
                let b = b.reopen()?;
                Ok(BigWigIntervalIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).convert_err()?),
                }
                .into_py(py))
            }
            #[cfg(feature = "remote")]
            BBIReadRaw::BigWigRemote(b) => {
                let b = b.reopen()?;
                Ok(BigWigIntervalIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).convert_err()?),
                }
                .into_py(py))
            }
            BBIReadRaw::BigWigFileLike(b) => {
                let b = b.reopen()?;
                Ok(BigWigIntervalIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).convert_err()?),
                }
                .into_py(py))
            }
            BBIReadRaw::BigBedFile(b) => {
                let b = b.reopen()?;
                Ok(BigBedEntriesIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).convert_err()?),
                }
                .into_py(py))
            }
            #[cfg(feature = "remote")]
            BBIReadRaw::BigBedRemote(b) => {
                let b = b.reopen()?;
                Ok(BigBedEntriesIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).convert_err()?),
                }
                .into_py(py))
            }
            BBIReadRaw::BigBedFileLike(b) => {
                let b = b.reopen()?;
                Ok(BigBedEntriesIterator {
                    iter: Box::new(b.get_interval_move(&chrom, start, end).convert_err()?),
                }
                .into_py(py))
            }
        }
    }

    /// Return the zoom records of a given range on a chromosome for a given
    /// zoom level.
    ///
    /// The result is an iterator of tuples. These tuples are in the format
    /// (start: int, end: int, summary: dict).
    ///
    /// Parameters
    /// ----------
    /// reduction_level : int
    ///     The zoom level to use, as a resolution in bases. Use the ``zooms``
    ///     method to get a list of available zoom levels.
    /// chrom : str
    ///     Name of the chromosome.
    /// start, end : int, optional
    ///     The range to get values for. If end is not provided, it defaults
    ///     to the length of the chromosome. If start is not provided, it
    ///     defaults to the beginning of the chromosome.
    ///
    /// Returns
    /// -------
    /// Iterator[tuple[int, int, dict]]
    ///     An iterator of tuples in the format (start: int, end: int,
    ///     summary: dict).
    ///
    /// Notes
    /// -----
    /// The summary dictionary contains the following keys
    ///
    /// - ``total_items``: The number of items in the interval.
    /// - ``bases_covered``: The number of bases covered by the interval.
    /// - ``min_val``: The minimum value in the interval.
    /// - ``max_val``: The maximum value in the interval.
    /// - ``sum``: The sum of all values in the interval.
    /// - ``sum_squares``: The sum of the squares of all values in the interval.
    ///
    /// For BigWigs, the summary statistics are derived from the unique
    /// **signal values** associated with each base in the interval.
    ///
    /// For BigBeds, the summary statistics instead are derived from the
    /// **number of BED intervals** overlapping each base in the interval.
    ///
    /// See Also
    /// --------
    /// zooms : Get a list of available zoom levels.
    /// records : Get the records of a given range on a chromosome.
    /// values : Get the values of a given range on a chromosome.
    #[pyo3(signature = (reduction_level, chrom, start=None, end=None))]
    fn zoom_records(
        &mut self,
        reduction_level: u32,
        chrom: String,
        start: Option<i32>,
        end: Option<i32>,
    ) -> PyResult<ZoomIntervalIterator> {
        let (start, end) = start_end(&self.bbi, &chrom, start, end)?;
        match &self.bbi {
            BBIReadRaw::Closed => return Err(BBIFileClosed::new_err("File is closed.")),
            BBIReadRaw::BigWigFile(b) => {
                let b = b.reopen()?;
                let iter = b
                    .get_zoom_interval_move(&chrom, start, end, reduction_level)
                    .convert_err()?;
                Ok(ZoomIntervalIterator {
                    iter: Box::new(iter),
                })
            }
            #[cfg(feature = "remote")]
            BBIReadRaw::BigWigRemote(b) => {
                let b = b.reopen()?;
                let iter = b
                    .get_zoom_interval_move(&chrom, start, end, reduction_level)
                    .convert_err()?;
                Ok(ZoomIntervalIterator {
                    iter: Box::new(iter),
                })
            }
            BBIReadRaw::BigWigFileLike(b) => {
                let b = b.reopen()?;
                let iter = b
                    .get_zoom_interval_move(&chrom, start, end, reduction_level)
                    .convert_err()?;
                Ok(ZoomIntervalIterator {
                    iter: Box::new(iter),
                })
            }
            BBIReadRaw::BigBedFile(b) => {
                let b = b.reopen()?;
                let iter = b
                    .get_zoom_interval_move(&chrom, start, end, reduction_level)
                    .convert_err()?;
                Ok(ZoomIntervalIterator {
                    iter: Box::new(iter),
                })
            }
            #[cfg(feature = "remote")]
            BBIReadRaw::BigBedRemote(b) => {
                let b = b.reopen()?;
                let iter = b
                    .get_zoom_interval_move(&chrom, start, end, reduction_level)
                    .convert_err()?;
                Ok(ZoomIntervalIterator {
                    iter: Box::new(iter),
                })
            }
            BBIReadRaw::BigBedFileLike(b) => {
                let b = b.reopen()?;
                let iter = b
                    .get_zoom_interval_move(&chrom, start, end, reduction_level)
                    .convert_err()?;
                Ok(ZoomIntervalIterator {
                    iter: Box::new(iter),
                })
            }
        }
    }

    /// Return the values of a given range on a chromosome as a numpy array.
    ///
    /// For BigWigs, the returned values or summary statistics are derived
    /// from the unique **signal values** associated with each base.
    ///
    /// For BigBeds, the returned values or summary statistics instead are
    /// derived from the **number of BED intervals** overlapping each base.
    ///
    /// Parameters
    /// ----------
    /// chrom : str
    ///     Name of the chromosome.  
    /// start, end : int, optional
    ///     The range to get values for. If end is not provided, it defaults
    ///     to the length of the chromosome. If start is not provided, it
    ///     defaults to the beginning of the chromosome.
    /// bins : int, optional
    ///     If provided, the query interval will be divided into equally spaced
    ///     bins and the values in each bin will be interpolated or summarized.
    ///     If not provided, the values will be returned for each base.
    /// summary : Literal["mean", "min", "max"], optional [default: "mean"]
    ///     The summary statistic to use. Currently supported statistics are
    ///     ``mean``, ``min``, and ``max``.
    /// exact : bool, optional [default: False]
    ///     If True and ``bins`` is specified, return exact summary statistic
    ///     values instead of interpolating from the optimal zoom level.
    ///     Default is False.
    /// missing : float, optional [default: 0.0]
    ///     Fill-in value for unreported data in valid regions. Default is 0.
    /// oob : float, optional [default: NaN]
    ///     Fill-in value for out-of-bounds regions. Default is NaN.
    /// arr : numpy.ndarray, optional
    ///     If provided, the values will be written to this array or array
    ///     view. The array must be of the correct size and type.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     The signal values of the bigwig or bigbed in the specified range.
    ///
    /// Notes
    /// -----
    /// A BigWig file encodes a step function, and the value at
    /// a base is given by the signal value of the unique interval that
    /// contains that base.
    ///
    /// A BigBed file encodes a collection of (possibly overlapping) intervals
    /// which may or may not be associated with quantitative scores. The
    /// "value" at given base used here summarizes the number of intervals
    /// overlapping that base, not any particular score.
    ///
    /// If a number of bins is requested and ``exact`` is False, the summarized
    /// data is interpolated from the closest available zoom level. If you
    /// need accurate summary data and are okay with small trade-off in speed,
    /// set ``exact`` to True.
    ///
    /// See Also
    /// --------
    /// records : Get the records of a given range on a chromosome.
    /// zoom_records : Get the zoom records of a given range on a chromosome.
    #[pyo3(
        signature = (chrom, start=None, end=None, bins=None, summary="mean".to_string(), exact=false, missing=0.0, oob=f64::NAN, arr=None),
        text_signature = r#"(chrom, start, end, bins=None, summary="mean", exact=False, missing=0.0, oob=..., arr=None)"#,
    )]
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
                return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
                    "Unrecognized summary. Only `mean`, `min`, and `max` are allowed."
                )));
            }
        };
        match &mut self.bbi {
            BBIReadRaw::Closed => return Err(BBIFileClosed::new_err("File is closed.")),
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

    /// Return the names of chromosomes in a BBI file and their lengths.  
    ///
    /// Parameters
    /// ----------
    /// chrom : str or None
    ///     The name of the chromosome to get the length of. If None, then a
    ///     dictionary of all chromosome sizes will be returned. If the
    ///     chromosome doesn't exist, returns None.
    ///
    /// Returns
    /// -------
    /// int or Dict[str, int] or None:
    ///     Chromosome length or a dictionary of chromosome lengths.
    #[pyo3(signature = (chrom=None))]
    fn chroms(&mut self, py: Python, chrom: Option<String>) -> PyResult<PyObject> {
        fn get_chrom_obj<B: bigtools::BBIRead>(
            b: &B,
            py: Python,
            chrom: Option<String>,
        ) -> PyResult<PyObject> {
            match chrom {
                Some(chrom) => {
                    let chrom_length = b
                        .chroms()
                        .into_iter()
                        .find(|c| c.name == chrom)
                        .ok_or_else(|| {
                            PyErr::new::<PyKeyError, _>(
                                "No chromosome found with the specified name",
                            )
                        })
                        .map(|c| c.length.to_object(py))?;
                    Ok(chrom_length)
                }
                None => {
                    let chrom_dict: PyObject = b
                        .chroms()
                        .into_iter()
                        .map(|c| (c.name.clone(), c.length))
                        .into_py_dict_bound(py)
                        .into();
                    Ok(chrom_dict)
                }
            }
        }

        match &self.bbi {
            BBIReadRaw::Closed => Err(BBIFileClosed::new_err("File is closed.")),
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

    /// Gets the average values from a bigWig over the entries of a bed file.
    ///
    /// Parameters
    /// ----------
    /// bed : str or Path
    ///     The path to the bed.
    /// names : bool or int, optional
    ///     If ``None``, then no name is returned and the return value is only the
    ///     statistics value (see the `stats` parameter).
    ///     
    ///     If ``True``, then each return value will be a 2-length tuple of the value
    ///     of column 4 and the statistics value.
    ///     
    ///     If ``False``, then each return value will be a 2-length tuple of the
    ///     interval in the format ``{chrom}:{start}-{end}`` and the statistics value.
    ///     
    ///     If ``0``, then each return value will match as if ``False`` was passed.
    ///   
    ///     If a ``1+``, then each return value will be a tuple of the value of
    ///     column of this parameter (1-based) and the statistics value.
    /// stats : str or List[str], optional
    ///     Calculate specific statistics for each bed entry.
    ///
    ///     If not specified, `mean` will be returned.
    ///
    ///     If ``"all"`` is specified, all summary statistics are returned in a named tuple.
    ///
    ///     If a single statistic is provided as a string, that statistic is returned
    ///     as a float or int depending on the statistic.
    ///
    ///     If a list of statistics are provided, a tuple is returned containing those
    ///     statistics, in order.
    ///
    ///     Possible statistics are:
    ///       - `size`: Size of bed entry (`int`)
    ///       - `bases`: Bases covered by bigWig (`int`)
    ///       - `sum`: Sum of values over all bases covered (`float`)
    ///       - `mean0`: Average over bases with non-covered bases counting as zeroes (`float`)
    ///       - `mean` or `None`: Average over just covered bases (`float`)
    ///       - `min`: Minimum over all bases covered (`float`)
    ///       - `max`: Maximum over all bases covered (`float`)
    ///
    /// Returns
    /// -------
    /// Generator of float or tuple.
    ///
    /// Notes
    /// -----
    /// If no ``name`` field is specified, returns a generator of statistics
    /// (either floats or tuples, as specified by the ``stats`` field).
    /// If a ``name`` column is specified, returns a generator of 2-length
    /// tuples of the form ``({name}, {average})``.
    /// Importantly, if the statistics value is itself a tuple, then that
    /// tuple will be **nested** as the second value of the outer tuple.
    #[pyo3(signature = (bed, names=None, stats=None))]
    fn average_over_bed(
        &mut self,
        py: Python,
        bed: Bound<'_, PyAny>,
        names: Option<PyObject>,
        stats: Option<PyObject>,
    ) -> PyResult<PyObject> {
        if self.is_bigbed() {
            // Minor thing: if current file is not a bigWig, don't even attempt to open the passed file
            return Err(PyErr::new::<exceptions::PyValueError, _>("Not a bigWig."));
        }
        let (name, usename) = {
            match names {
                Some(names) => match names.extract::<bool>(py) {
                    Ok(true) => (Name::Column(3), true),
                    Ok(false) => (Name::Interval, true),
                    Err(_) => match names.extract::<isize>(py) {
                        Ok(col) => match col {
                            0 => (Name::None, true),
                            1.. => (Name::Column((col - 1) as usize), true),
                            _ => {
                                return Err(PyErr::new::<exceptions::PyValueError, _>(
                                    "Invalid names argument. Must be >= 0.",
                                ));
                            }
                        },
                        Err(_) => {
                            return Err(PyErr::new::<exceptions::PyValueError, _>("Invalid names argument. Should be either `None`, a `bool`, or an `int`"));
                        }
                    },
                },
                None => (Name::Interval, false),
            }
        };
        let stats = {
            match stats {
                Some(stats) => match stats.downcast_bound::<PyString>(py) {
                    Ok(stat) => {
                        let stat = stat.to_str()?;
                        if stat.eq_ignore_ascii_case("all") {
                            None
                        } else {
                            let stat = match BigWigAverageOverBedStatistics::from_str(stat) {
                                Some(stat) => stat,
                                None => {
                                    return Err(PyErr::new::<exceptions::PyValueError, _>("Invalid type argument. Should be either `None` or \"size\", \"bases\", \"sum\", \"mean0\", \"min\", \"max\" or \"all\""));
                                }
                            };
                            Some(vec![stat])
                        }
                    }
                    Err(_) => match stats.downcast_bound::<PyList>(py) {
                        Ok(stats) => {
                            let mut ret_stats = Vec::with_capacity(stats.len());
                            for stat in stats.into_iter() {
                                let stat = stat.downcast::<PyString>();
                                let Ok(stat) = stat else {
                                    return Err(PyErr::new::<exceptions::PyValueError, _>("Invalid type argument. Should be either `None` or \"all\" or a stat or list of stats."));
                                };
                                let stat = match BigWigAverageOverBedStatistics::from_str(
                                    stat.to_str()?,
                                ) {
                                    Some(stat) => stat,
                                    None => {
                                        return Err(PyErr::new::<exceptions::PyValueError, _>("Invalid type argument. Should be either `None` or \"size\", \"bases\", \"sum\", \"mean0\", \"min\", \"max\" or \"all\""));
                                    }
                                };
                                ret_stats.push(stat);
                            }
                            Some(ret_stats)
                        }
                        Err(_) => {
                            return Err(PyErr::new::<exceptions::PyValueError, _>("Invalid type argument. Should be either `None` or \"all\" or a stat or list of stats."));
                        }
                    },
                },
                None => Some(vec![BigWigAverageOverBedStatistics::Mean]),
            }
        };
        let bed = if let Ok(bed) = bed.downcast::<PyString>() {
            bed.clone()
        } else {
            let path_class = py.import_bound("pathlib")?.getattr("Path")?;
            // If pathlib.Path, convert to string and try to open
            if bed.is_instance(&path_class)? {
                bed.str()?
            } else {
                return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
                    "Unknown argument for `path`. Not a string or Path object.",
                )));
            }
        };
        let bedin = BufReader::new(File::open(bed.to_str()?)?);

        let module = PyModule::import_bound(py, "pybigtools")?;
        let summary_statistics = module.getattr("SummaryStatistics")?.to_object(py);
        let res = match &mut self.bbi {
            BBIReadRaw::Closed => return Err(BBIFileClosed::new_err("File is closed.")),
            BBIReadRaw::BigWigFile(b) => {
                let b = b.reopen()?;
                let iter = Box::new(bigwig_average_over_bed(bedin, b, name));
                BigWigAverageOverBedEntriesIterator {
                    iter,
                    usename,
                    stats,
                    summary_statistics,
                }
                .into_py(py)
            }
            #[cfg(feature = "remote")]
            BBIReadRaw::BigWigRemote(b) => {
                let b = b.reopen()?;
                let iter = Box::new(bigwig_average_over_bed(bedin, b, name));
                BigWigAverageOverBedEntriesIterator {
                    iter,
                    usename,
                    stats,
                    summary_statistics,
                }
                .into_py(py)
            }
            BBIReadRaw::BigWigFileLike(b) => {
                let b = b.reopen()?;
                let iter = Box::new(bigwig_average_over_bed(bedin, b, name));
                BigWigAverageOverBedEntriesIterator {
                    iter,
                    usename,
                    stats,
                    summary_statistics,
                }
                .into_py(py)
            }
            BBIReadRaw::BigBedFile(_) | BBIReadRaw::BigBedFileLike(_) => {
                return Err(PyErr::new::<exceptions::PyValueError, _>("Not a bigWig."));
            }
            #[cfg(feature = "remote")]
            BBIReadRaw::BigBedRemote(_) => {
                return Err(PyErr::new::<exceptions::PyValueError, _>("Not a bigWig."))
            }
        };

        Ok(res)
    }

    fn close(&mut self) {
        self.bbi = BBIReadRaw::Closed;
    }

    fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    fn __exit__(
        slf: Py<Self>,
        py: Python<'_>,
        _exc_type: PyObject,
        _exc_value: PyObject,
        _exc_traceback: PyObject,
    ) {
        slf.borrow_mut(py).bbi = BBIReadRaw::Closed;
    }

    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        match &self.bbi {
            BBIReadRaw::Closed | BBIReadRaw::BigWigFile(_) | BBIReadRaw::BigBedFile(_) => {}
            #[cfg(feature = "remote")]
            BBIReadRaw::BigWigRemote(_) | BBIReadRaw::BigBedRemote(_) => {}
            BBIReadRaw::BigWigFileLike(b) => {
                visit.call(b.inner_read().inner_read().inner.deref())?
            }
            BBIReadRaw::BigBedFileLike(b) => {
                visit.call(b.inner_read().inner_read().inner.deref())?
            }
        }
        Ok(())
    }

    fn __clear__(&mut self) {
        self.close()
    }
}

#[pyclass(module = "pybigtools")]
struct ZoomIntervalIterator {
    iter: Box<dyn Iterator<Item = Result<ZoomRecord, _BBIReadError>> + Send>,
}

#[pymethods]
impl ZoomIntervalIterator {
    fn __iter__(slf: PyRefMut<Self>) -> PyResult<Py<ZoomIntervalIterator>> {
        Ok(slf.into())
    }

    fn __next__(mut slf: PyRefMut<Self>) -> PyResult<Option<(u32, u32, PyObject)>> {
        slf.iter
            .next()
            .transpose()
            .map(|o| {
                o.map(|v| {
                    let summary = [
                        ("total_items", v.summary.total_items.to_object(slf.py())),
                        ("bases_covered", v.summary.bases_covered.to_object(slf.py())),
                        ("min_val", v.summary.min_val.to_object(slf.py())),
                        ("max_val", v.summary.max_val.to_object(slf.py())),
                        ("sum", v.summary.sum.to_object(slf.py())),
                        ("sum_squares", v.summary.sum_squares.to_object(slf.py())),
                    ]
                    .into_py_dict_bound(slf.py())
                    .to_object(slf.py());
                    (v.start, v.end, summary)
                })
            })
            .convert_err()
    }
}

/// An iterator for intervals in a bigWig.  
///
/// It returns only values that exist in the bigWig, skipping any missing
/// intervals.
#[pyclass(module = "pybigtools")]
pub struct BigWigIntervalIterator {
    iter: Box<dyn Iterator<Item = Result<Value, _BBIReadError>> + Send>,
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

/// An iterator for the entries in a bigBed.
#[pyclass(module = "pybigtools")]
pub struct BigBedEntriesIterator {
    iter: Box<dyn Iterator<Item = Result<BedEntry, _BBIReadError>> + Send>,
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
            PyTuple::new_bound::<PyObject, _>(py, elements.into_iter()).to_object(py),
        ))
    }
}
