# ruff: noqa: F403
# ruff: noqa: F405
# This is here to make sure that source installs with pip work correctly
import warnings as _warnings

from pybigtools.pybigtools import *
from pybigtools import pybigtools as bt

__doc__ = bt.__doc__
if hasattr(bt, "__all__"):
    __all__ = bt.__all__

_RustBBIReader = bt.BBIReader
_rust_open = bt.open


# Sentinel distinguishing "user did not pass" from "user passed None" in `values`.
_UNSET = object()


class BBIReader:
    """Python-side wrapper around the Rust ``BBIReader`` that stages a
    deprecation cycle on the ``values()`` parameters. All other attributes
    are delegated to the underlying Rust object.
    """

    __slots__ = ("_rust",)

    def __init__(self, rust_reader):
        object.__setattr__(self, "_rust", rust_reader)

    def __getattr__(self, name):
        return getattr(self._rust, name)

    def __enter__(self):
        self._rust.__enter__()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return self._rust.__exit__(exc_type, exc_value, traceback)

    def values(
        self,
        chrom,
        start=None,
        end=None,
        bins=None,
        summary="mean",
        exact=False,
        uncovered=None,
        oob=float("nan"),
        fillna=_UNSET,
        missing=_UNSET,
        arr=None,
    ):
        if missing is not _UNSET:
            _warnings.warn(
                "`missing` is deprecated and will be removed in a future "
                "release; use `fillna` instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if fillna is _UNSET:
                fillna = missing
        if fillna is _UNSET:
            _warnings.warn(
                "The default behavior of `values()` has changed: empty bins / "
                "uncovered positions are now returned as NaN instead of "
                "filled with 0. Pass `fillna=0` to keep the previous "
                "behavior, or `fillna=None` to silence this warning. See "
                "also `uncovered` to include uncovered bases in summary "
                "statistic calculations.",
                DeprecationWarning,
                stacklevel=2,
            )
            fillna = None
        return self._rust.values(
            chrom,
            start,
            end,
            bins=bins,
            summary=summary,
            exact=exact,
            uncovered=uncovered,
            oob=oob,
            fillna=fillna,
            arr=arr,
        )


def open(path, mode="r"):  # noqa: A001 - intentional shadow of builtin
    return BBIReader(_rust_open(path, mode))


del bt
