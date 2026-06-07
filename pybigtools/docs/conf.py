"""Sphinx configuration for the pybigtools documentation.

The API reference is built with autodoc by importing the compiled extension,
so the build needs pybigtools installed (see .readthedocs.yaml). Docstrings are
NumPy-style and rendered via napoleon; signatures come from the extension's
``__text_signature__`` metadata.
"""

import datetime
import importlib.metadata

import pybigtools

project = "pybigtools"
author = "Bigtools contributors"
copyright = f"2024-{datetime.date.today().year}, {author}"

try:
    release = importlib.metadata.version("pybigtools")
except importlib.metadata.PackageNotFoundError:
    release = pybigtools.__version__
version = release

pybigtools_version = pybigtools.__version__
bigtools_version = pybigtools.__core_version__

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
html_title = "Bigtools"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
templates_path = ["_templates"]
html_theme_options = {
    "light_logo": "bigtools-logo.svg",
    "dark_logo": "bigtools-logo-dark.svg",
}
# Exposed to the sidebar brand template.
html_context = {
    "pybigtools_version": pybigtools_version,
    "bigtools_version": bigtools_version,
}

# -- MyST --------------------------------------------------------------------
myst_enable_extensions = ["colon_fence", "deflist"]

exclude_patterns = ["_build"]
