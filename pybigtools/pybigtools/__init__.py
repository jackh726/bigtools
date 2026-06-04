# ruff: noqa: F403
# ruff: noqa: F405
# This is here to make sure that source installs with pip work correctly
import inspect as _inspect
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
    """Interface for reading a BigWig or BigBed file.

    Returned by :func:`open` in read mode. Use the methods below to query
    chromosomes, intervals, records, zoom levels, and summary statistics, or
    extract values as NumPy arrays. Supports the context-manager protocol.
    """

    # Implementation note: this is a thin Python wrapper around the Rust
    # ``BBIReader`` that stages a deprecation cycle on the ``values()``
    # parameters; all other attributes are delegated to the Rust object.

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


# ---------------------------------------------------------------------------
# Make the public wrapper objects fully introspectable.
#
# ``BBIReader`` delegates most calls to the underlying Rust reader via
# ``__getattr__``, which means those methods are invisible to ``help()``, IDEs,
# and documentation tools. Here we surface them as real methods/properties that
# forward to the Rust object, carrying the upstream docstrings and signatures.
# This is purely cosmetic/introspective; behavior is unchanged.
# ---------------------------------------------------------------------------


def _signature_from_text(text_signature, *, prepend_self=False):
    """Build an ``inspect.Signature`` from a C ``__text_signature__`` string."""
    source = f"def _f{text_signature}: ..."
    namespace = {}
    exec(source, {"__builtins__": __builtins__}, namespace)  # noqa: S102
    sig = _inspect.signature(namespace["_f"])
    if prepend_self:
        params = [
            _inspect.Parameter("self", _inspect.Parameter.POSITIONAL_OR_KEYWORD),
            *sig.parameters.values(),
        ]
        sig = sig.replace(parameters=params)
    return sig


def _make_delegator(name, rust_method):
    def method(self, *args, **kwargs):
        return getattr(self._rust, name)(*args, **kwargs)

    method.__name__ = name
    method.__qualname__ = f"BBIReader.{name}"
    method.__doc__ = rust_method.__doc__
    try:
        method.__signature__ = _inspect.signature(rust_method)
    except (TypeError, ValueError):
        text = getattr(rust_method, "__text_signature__", None)
        if text:
            method.__signature__ = _signature_from_text(text, prepend_self=True)
    return method


for _name in (
    "chroms",
    "info",
    "zooms",
    "records",
    "zoom_records",
    "average_over_bed",
    "sql",
    "close",
):
    setattr(BBIReader, _name, _make_delegator(_name, getattr(_RustBBIReader, _name)))

for _name in ("is_bigwig", "is_bigbed"):
    _rust_attr = getattr(_RustBBIReader, _name)
    setattr(
        BBIReader,
        _name,
        property(
            (lambda n: lambda self: getattr(self._rust, n))(_name),
            doc=_rust_attr.__doc__,
        ),
    )

# The deprecation wrapper above redefines ``values``; give it the upstream
# docstring and the clean Rust signature for display purposes.
BBIReader.values.__doc__ = _RustBBIReader.values.__doc__
BBIReader.values.__signature__ = _signature_from_text(
    _RustBBIReader.values.__text_signature__, prepend_self=True
)

open.__doc__ = _rust_open.__doc__

del bt
