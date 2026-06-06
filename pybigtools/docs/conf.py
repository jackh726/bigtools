"""Sphinx configuration for the pybigtools documentation.

The API reference is built with autodoc by importing the compiled extension,
so the build needs pybigtools installed (see .readthedocs.yaml). Docstrings are
NumPy-style and rendered via napoleon; signatures come from the extension's
``__text_signature__`` metadata.
"""

import importlib.metadata

project = "pybigtools"
author = "Bigtools contributors"
copyright = "2024-2026, " + author

try:
    release = importlib.metadata.version("pybigtools")
except importlib.metadata.PackageNotFoundError:
    release = ""
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "myst_parser",
]

# -- Autodoc / napoleon -----------------------------------------------------
autodoc_member_order = "groupwise"
add_module_names = False
python_use_unqualified_type_names = True

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True

# -- Cross-project references ------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}

# -- HTML output -------------------------------------------------------------
html_theme = "furo"
html_title = f"pybigtools {release}".strip()

# -- MyST --------------------------------------------------------------------
myst_enable_extensions = ["colon_fence", "deflist"]

exclude_patterns = ["_build"]
