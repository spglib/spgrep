# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "spgrep"
copyright = "2022, Kohei Shinohara"
author = "Kohei Shinohara"

# https://github.com/pypa/setuptools_scm/
from importlib.metadata import version

release = version("spgrep")
# for example take major/minor
version = ".".join(release.split(".")[:3])

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.mermaid",
    "nbsphinx",
    "myst_parser",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["README.md"]

# The suffix(es) of source filenames.
source_suffix = [".rst"]

autoclass_content = "both"
autodoc_typehints = "signature"
autodoc_member_order = "bysource"
autodoc_type_aliases = {}

# napoleon_type_aliases = {}
napoleon_use_rtype = True
napoleon_use_ivar = True

# https://pypi.org/project/sphinxcontrib-bibtex/
bibtex_bibfiles = ["references.bib"]

# MyST
myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "tasklist",
]
myst_dmath_double_inline = True

# nbsphinx and sphinxcontrib.mermaid are conflicted.
# So, we need to use mermaid CLI instead of a raw HTML output.
# https://github.com/mgaitan/sphinxcontrib-mermaid/issues/74
mermaid_output_format = "svg"
mermaid_verbose = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
html_title = project + " " + version
html_theme_options = {
    "navigation_with_keys": True,
}

# hide sphinx footer
html_show_sphinx = False
html_show_sourcelink = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
