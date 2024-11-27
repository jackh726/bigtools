# This is here to make sure that source installs with pip work correctly

from .pybigtools import *

__doc__ = pybigtools.__doc__
if hasattr(pybigtools, "__all__"):
    __all__ = pybigtools.__all__
