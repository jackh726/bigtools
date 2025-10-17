#![allow(non_snake_case)]

use pyo3::prelude::*;
use pyo3::{create_exception, exceptions};

create_exception!(
    pybigtools,
    BBIFileClosed,
    exceptions::PyException,
    "BBI File is closed."
);
create_exception!(
    pybigtools,
    BBIReadError,
    exceptions::PyException,
    "Error reading BBI file."
);

trait ToPyErr {
    fn to_py_err(self) -> PyErr;
}

impl ToPyErr for bigtools::BBIReadError {
    fn to_py_err(self) -> PyErr {
        PyErr::new::<BBIReadError, _>(format!("{}", self))
    }
}
impl ToPyErr for bigtools::ZoomIntervalError {
    fn to_py_err(self) -> PyErr {
        match self {
            bigtools::ZoomIntervalError::ReductionLevelNotFound => {
                PyErr::new::<exceptions::PyKeyError, _>(format!(
                    "The passed reduction level was not found"
                ))
            }
            _ => PyErr::new::<BBIReadError, _>(format!("{}", self)),
        }
    }
}

pub trait ConvertResult<T> {
    fn convert_err(self) -> Result<T, PyErr>;
}
impl<T, E: ToPyErr> ConvertResult<T> for Result<T, E> {
    fn convert_err(self) -> Result<T, PyErr> {
        self.map_err(|e| e.to_py_err())
    }
}
