use std::io;

use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BedEntry, BigBedWrite as BigBedWriteRaw, BigWigWrite as BigWigWriteRaw, Value};
use pyo3::exceptions::{self, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyFloat, PyInt, PyIterator, PyString, PyTuple};
use tokio::runtime;

use crate::errors::BBIFileClosed;

/// Interface for writing to a BigWig file.
#[pyclass(module = "pybigtools")]
pub struct BigWigWrite {
    pub bigwig: Option<String>,
}

#[pymethods]
impl BigWigWrite {
    /// Write values to the BigWig file.
    ///
    /// The underlying file will be closed automatically when the function
    /// completes (and no other operations will be able to be performed).
    ///
    /// Parameters
    /// ----------
    /// chroms : Dict[str, int]
    ///     A dictionary with keys as chromosome names and values as their
    ///     length.
    /// vals : Iterable[tuple[str, int, int, float]]
    ///     An iterable with values that represents each value to write in the
    ///     format (chromosome, start, end, value).
    ///
    /// Notes
    /// -----
    /// The underlying file will be closed automatically when the function
    /// completes, and no other operations will be able to be performed.  
    fn write(&mut self, py: Python, chroms: &Bound<'_, PyDict>, vals: Py<PyAny>) -> PyResult<()> {
        let runtime = runtime::Builder::new_multi_thread()
            .worker_threads(
                std::thread::available_parallelism()
                    .map(|c| c.into())
                    .unwrap_or(1),
            )
            .build()
            .expect("Unable to create thread pool.");

        let chrom_map = chroms
            .into_iter()
            .map(|(key, val)| {
                let chrom: String = key.downcast::<PyString>()?.to_str().unwrap().to_owned();
                let length: u32 = val.downcast::<PyInt>()?.to_object(py).extract(py).unwrap();
                Ok((chrom, length))
            })
            .collect::<Result<std::collections::HashMap<String, u32>, pyo3::PyErr>>()?;

        let bigwig = self
            .bigwig
            .take()
            .ok_or_else(|| PyErr::new::<BBIFileClosed, _>("Can only write once."))?;
        let bigwig = BigWigWriteRaw::create_file(bigwig, chrom_map).map_err(|e| {
            PyErr::new::<exceptions::PyException, _>(format!(
                "Error occured when creating file: {}",
                e
            ))
        })?;

        struct IterError(String);
        struct Iter {
            inner: PyObject,
        }
        impl Iterator for Iter {
            type Item = Result<(String, Value), IterError>;
            fn next(&mut self) -> Option<Self::Item> {
                // We have to reacquire the gil for each iteration
                Python::with_gil(|py| {
                    let mut iter: Bound<'_, PyIterator> = match self.inner.downcast_bound(py) {
                        Ok(o) => o.clone(),
                        Err(_) => {
                            return Some(Err(IterError(format!(
                                "Passed value for `val` is not iterable."
                            ))))
                        }
                    };
                    let next: Result<(String, Value), pyo3::PyErr> = match iter.next()? {
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
                let inner_obj = vals.bind(py);
                match PyIterator::from_bound_object(inner_obj) {
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
            let data = BedParserStreamingIterator::wrap_iter(vals_iter_raw, true);
            match bigwig.write(data, runtime) {
                Err(e) => println!("{}", e),
                Ok(_) => {}
            }
            Ok(())
        })
    }

    /// Close the file.
    ///
    /// No other operations will be allowed after it is closed. This is done
    /// automatically after write is performed.
    fn close(&mut self) -> PyResult<()> {
        self.bigwig.take();
        Ok(())
    }
}

/// Interface for writing to a BigBed file.
#[pyclass(module = "pybigtools")]
pub struct BigBedWrite {
    pub bigbed: Option<String>,
}

#[pymethods]
impl BigBedWrite {
    /// Write values to the BigBed file.
    ///
    /// The underlying file will be closed automatically when the function
    /// completes (and no other operations will be able to be performed).
    ///
    /// Parameters
    /// ----------
    /// chroms : Dict[str, int]
    ///     A dictionary with keys as chromosome names and values as their
    ///     length.
    /// vals : Iterable[tuple[str, int, int, str]]
    ///     An iterable with values that represents each value to write in the
    ///     format (chromosome, start, end, rest). The ``rest`` string should
    ///     consist of tab-delimited fields.  
    ///
    /// Notes
    /// -----
    /// The underlying file will be closed automatically when the function
    /// completes, and no other operations will be able to be performed.
    #[pyo3(signature = (chroms, vals, autosql=None))]
    fn write(
        &mut self,
        py: Python,
        chroms: &Bound<'_, PyDict>,
        vals: Py<PyAny>,
        autosql: Option<Bound<'_, PyString>>,
    ) -> PyResult<()> {
        let runtime = runtime::Builder::new_multi_thread()
            .worker_threads(
                std::thread::available_parallelism()
                    .map(|c| c.into())
                    .unwrap_or(1),
            )
            .build()
            .expect("Unable to create thread pool.");

        let chrom_map = chroms
            .into_iter()
            .map(|(key, val)| {
                let chrom: String = key.downcast::<PyString>()?.to_str().unwrap().to_owned();
                let length: u32 = val.downcast::<PyInt>()?.to_object(py).extract(py).unwrap();
                Ok((chrom, length))
            })
            .collect::<Result<std::collections::HashMap<String, u32>, pyo3::PyErr>>()?;

        let bigbed = self
            .bigbed
            .take()
            .ok_or_else(|| PyErr::new::<BBIFileClosed, _>("File already closed."))?;
        let mut bigbed = BigBedWriteRaw::create_file(bigbed, chrom_map).map_err(|e| {
            PyErr::new::<exceptions::PyException, _>(format!(
                "Error occured when creating file: {}",
                e
            ))
        })?;
        if let Some(autosql) = autosql {
            bigbed.autosql = Some(autosql.str()?.to_string());
        }

        struct IterError(String);
        struct Iter {
            inner: PyObject,
        }
        impl Iterator for Iter {
            type Item = Result<(String, BedEntry), IterError>;
            fn next(&mut self) -> Option<Self::Item> {
                // We have to reacquire the gil for each iteration
                Python::with_gil(|py| {
                    let mut iter: Bound<'_, PyIterator> = match self.inner.downcast_bound(py) {
                        Ok(o) => o.clone(),
                        Err(_) => {
                            return Some(Err(IterError(format!(
                                "Passed value for `val` is not iterable."
                            ))))
                        }
                    };
                    let next: Result<(String, BedEntry), pyo3::PyErr> = match iter.next()? {
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
                                    .get_item(3)
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
                let inner_obj = vals.bind(py);
                match PyIterator::from_bound_object(inner_obj) {
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
            let data = BedParserStreamingIterator::wrap_iter(vals_iter_raw, true);
            match bigbed.write(data, runtime) {
                Err(e) => {
                    println!("{}", e)
                }
                Ok(_) => {}
            }
            Ok(())
        })
    }

    /// Close the file.
    ///
    /// No other operations will be allowed after it is closed. This is done
    /// automatically after write is performed.
    fn close(&mut self) -> PyResult<()> {
        self.bigbed.take();
        Ok(())
    }
}
