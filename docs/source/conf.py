# -*- coding: utf-8 -*-
import sys
from pathlib import Path

import cellrank  # noqa

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent.parent))

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
needs_sphinx = "2.0"

# -- Project information -----------------------------------------------------

project = "CellRank"
copyright = "2019, Marius Lange, Michal Klein, Juan Luis Restrepo Lopez"
author = "Marius Lange, Michal Klein, Juan Luis Restrepo Lopez"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx.ext.intersphinx",
    "sphinx_paramlinks",
    "sphinx.ext.autosummary",
]

intersphinx_mapping = dict(
    anndata=("https://anndata.readthedocs.io/en/stable/", None),
    scanpy=("https://scanpy.readthedocs.io/en/stable/", None),
    scvelo=("https://scvelo.readthedocs.io/", None),
    python=("https://docs.python.org/3", None),
    numpy=("https://docs.scipy.org/doc/numpy/", None),
    scipy=("https://docs.scipy.org/doc/scipy/reference/", None),
    networkx=("https://networkx.github.io/documentation/stable/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    statsmodels=("https://www.statsmodels.org/stable/", None),
    matplotlib=("https://matplotlib.org/", None),
    joblib=("https://joblib.readthedocs.io/en/latest/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
    seaborn=("https://seaborn.pydata.org/", None),
    # TODO: add msmtools once the docs are up
)

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"
pygments_style = "sphinx"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

autosummary_generate = True
autodoc_member_order = "bysource"
# autodoc_default_flags = ['members']
typehints_fully_qualified = False
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = dict(navigation_depth=4, logo_only=True)
html_context = dict(
    display_github=True,  # Integrate GitHub
    github_user="theislab",  # Username
    github_repo="cellrank",  # Repo name
    github_version="master",  # Version
    conf_py_path="/docs/source/",
)  # Path in the checkout to the docs root
html_show_sphinx = False


def setup(app):
    app.add_stylesheet("css/custom.css")
