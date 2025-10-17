use std::io;
use std::path::Path;

#[cfg(feature = "remote")]
use bigtools::utils::file::remote_file::RemoteFile;
use bigtools::utils::reopen::Reopen;
use bigtools::{BigBedRead as BigBedReadRaw, BigWigRead as BigWigReadRaw, GenericBBIRead};
use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyString};
use pyo3::wrap_pyfunction;
use url::Url;

mod errors;
mod file_like;
mod reader;
mod utils;
mod writer;

use errors::{BBIFileClosed, BBIReadError};
use file_like::PyFileLikeObject;
use reader::{BBIRead, BigBedEntriesIterator, BigWigIntervalIterator};
use utils::BBIReadRaw;
use writer::{BigBedWrite, BigWigWrite};

impl Reopen for PyFileLikeObject {
    fn reopen(&self) -> io::Result<Self> {
        Ok(self.clone())
    }
}

/// Open a BigWig or BigBed file for reading or writing.
///
/// Parameters
/// ----------
/// path_url_or_file_like : str or file-like object
///     The path to a file or an http url for a remote file as a string, or
///     a Python file-like object with ``read`` and ``seek`` methods.
/// mode : Literal["r", "w"], optional [default: "r"]
///     The mode to open the file in. If not provided, it will default to read.
///     "r" will open a bigWig/bigBed for reading but will not allow writing.
///     "w" will open a bigWig/bigBed for writing but will not allow reading.
///
/// Returns
/// -------
/// BigWigWrite or BigBedWrite or BBIRead
///     The object for reading or writing the BigWig or BigBed file.
///
/// Notes
/// -----
/// For writing, only a file path is currently accepted.
///
/// If passing a file-like object, concurrent reading of different intervals
/// is not supported and may result in incorrect behavior.
#[pyfunction]
#[pyo3(signature = (path_url_or_file_like, mode=None))]
fn open(
    py: Python,
    path_url_or_file_like: &Bound<'_, PyAny>,
    mode: Option<String>,
) -> PyResult<PyObject> {
    let iswrite = match &mode {
        Some(mode) if mode == "w" => true,
        Some(mode) if mode == "r" => false,
        None => false,
        Some(mode) => {
            return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
                "Invalid mode: `{}`",
                mode
            )));
        }
    };

    // If string, might be path or url like
    if let Ok(string_ref) = path_url_or_file_like.downcast::<PyString>() {
        return open_path_or_url(py, string_ref.to_str().unwrap().to_owned(), iswrite);
    }

    // If pathlib.Path, convert to string and try to open
    let path_class = py.import_bound("pathlib")?.getattr("Path")?;
    if path_url_or_file_like.is_instance(&path_class)? {
        let path_string = path_url_or_file_like.str()?.to_string();
        return open_path_or_url(py, path_string, iswrite);
    }

    if iswrite {
        return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
            "Writing only supports file names",
        )));
    }
    let file_like = match PyFileLikeObject::new(path_url_or_file_like, true, false, true) {
        Ok(file_like) => file_like,
        Err(_) => return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
            "Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.",
        ))),
    };
    let read = match GenericBBIRead::open(file_like) {
        Ok(GenericBBIRead::BigWig(bigwig)) => BBIRead {
            bbi: BBIReadRaw::BigWigFileLike(bigwig.cached()),
        }
        .into_py(py),
        Ok(GenericBBIRead::BigBed(bigbed)) => BBIRead {
            bbi: BBIReadRaw::BigBedFileLike(bigbed.cached()),
        }
        .into_py(py),
        Err(e) => {
            return Err(PyErr::new::<BBIReadError, _>(format!(
            "File-like object is not a bigWig or bigBed. Or there was just a problem reading: {e}",
        )))
        }
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
            return Err(PyErr::new::<exceptions::PyValueError, _>(format!("Invalid file type. Must be either a bigWig (.bigWig, .bw) or bigBed (.bigBed, .bb).")));
        }
    };
    let res = if iswrite {
        match Url::parse(&path_url_or_file_like) {
            Ok(_) => {
                return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
                    "Invalid file path. Writing does not support urls."
                )))
            }
            Err(_) => {}
        }
        match extension.as_ref() {
            "bw" | "bigWig" | "bigwig" => BigWigWrite {
                bigwig: Some(path_url_or_file_like),
            }
            .into_py(py),
            "bb" | "bigBed" | "bigbed" => BigBedWrite {
                bigbed: Some(path_url_or_file_like),
            }
            .into_py(py),
            _ => {
                return Err(PyErr::new::<exceptions::PyValueError, _>(format!("Invalid file type. Must be either a bigWig (.bigWig, .bw) or bigBed (.bigBed, .bb).")));
            }
        }
    } else {
        let isfile = Path::new(&path_url_or_file_like).exists();
        if !isfile {
            match Url::parse(&path_url_or_file_like) {
                Ok(_) => {}
                Err(_) => {
                    return Err(PyErr::new::<exceptions::PyValueError, _>(format!(
                        "Invalid file path. The file does not exist and it is not a url."
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
                            return Err(PyErr::new::<BBIReadError, _>(format!(
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
                            return Err(PyErr::new::<BBIReadError, _>(format!(
                                "Error opening bigWig."
                            )))
                        }
                    }

                    #[cfg(not(feature = "remote"))]
                    {
                        return Err(PyErr::new::<exceptions::PyOSError, _>(format!(
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
                            return Err(PyErr::new::<BBIReadError, _>(format!(
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
                            return Err(PyErr::new::<BBIReadError, _>(format!(
                                "Error opening bigBed."
                            )))
                        }
                    }

                    #[cfg(not(feature = "remote"))]
                    {
                        return Err(PyErr::new::<exceptions::PyOSError, _>(format!(
                            "Builtin support for remote files is not supported on this platform."
                        )));
                    }
                }
            }
            _ => {
                return Err(PyErr::new::<exceptions::PyValueError, _>(format!("Invalid file type. Must be either a bigWig (.bigWig, .bw) or bigBed (.bigBed, .bb).")));
            }
        }
    };
    Ok(res)
}

/// Read and write Big Binary Indexed (BBI) file types: BigWig and BigBed.
#[pymodule]
fn pybigtools(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    m.add_wrapped(wrap_pyfunction!(open))?;

    m.add_class::<BBIRead>()?;
    m.add_class::<BigWigWrite>()?;
    m.add_class::<BigBedWrite>()?;
    m.add_class::<BigWigIntervalIterator>()?;
    m.add_class::<BigBedEntriesIterator>()?;

    m.add("BBIFileClosed", m.py().get_type_bound::<BBIFileClosed>())?;
    m.add("BBIReadError", m.py().get_type_bound::<BBIReadError>())?;

    let collections = PyModule::import_bound(py, "collections")?;
    let namedtuple = collections.getattr("namedtuple")?;

    let summary_statistics = namedtuple.call(
        (
            "SummaryStatistics".to_object(py),
            ("size", "bases", "sum", "mean0", "mean", "min", "max").to_object(py),
        ),
        None,
    )?;
    m.add("SummaryStatistics", summary_statistics)?;

    Ok(())
}
