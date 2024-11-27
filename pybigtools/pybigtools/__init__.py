# ruff: noqa: F403
# ruff: noqa: F405
# This is here to make sure that source installs with pip work correctly
from pybigtools.pybigtools import *
from pybigtools import pybigtools as bt

__doc__ = bt.__doc__
if hasattr(bt, "__all__"):
    __all__ = bt.__all__
del bt
