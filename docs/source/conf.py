# -*- coding: utf-8 -*-
import os
import sys
import logging
import subprocess
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import urlretrieve

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, os.path.abspath("_ext"))

# this must be called prior to importing CellRank
if not os.path.exists(os.path.join(sys.path[1], "cellrank", "_vendor")):
    config_path = os.path.join(HERE.parent.parent, "vendorize.toml")
    print(f"Running vendorize using config: {config_path}")
    subprocess.run(["python-vendorize", config_path])


import cellrank  # noqa NOQA


logger = logging.getLogger(__name__)

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
needs_sphinx = "3.0"

notebooks_url = "https://github.com/theislab/cellrank_notebooks/raw/master/tutorials/"
for nb in ["pancreas_basic.ipynb", "pancreas_advanced.ipynb"]:
    try:
        url = urljoin(notebooks_url, nb)
        urlretrieve(url, nb)
    except Exception as e:
        logger.error(f"Unable to retrieve notebook: `{url}`. Reason: `{e}`.")


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
    "nbsphinx",
    "sphinx_copybutton",
    "sphinx_last_updated_by_git",
    "edit_on_github",
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
source_suffix = [".rst", ".ipynb"]
master_doc = "index"
pygments_style = "sphinx"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["**.ipynb_checkpoints"]


# -- Notebooks
nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'png', 'pdf'}",  # correct figure resize
    "--InlineBackend.rc={'figure.dpi': 96}",
]

nbsphinx_prolog = r"""
{% set docname = 'tutorials/' + env.doc2path(env.docname, base=None) %}
.. raw:: html

    <div class="note">
      Interactive version
      <a href="https://mybinder.org/v2/gh/theislab/cellrank_notebooks/{{ env.config.release|e }}?filepath={{ docname|e }}"><img alt="Binder badge" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a>
    </div>
"""

release = "master"


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

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
github_repo = "cellrank"  # sets the html_context
github_nb_repo = "cellrank_notebooks"
html_show_sphinx = False


def setup(app):
    app.add_stylesheet("css/custom.css")
