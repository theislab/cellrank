# based on: https://github.com/scverse/scvi-tools/blob/master/docs/conf.py
import os
import sys
import logging
from pathlib import Path
from datetime import datetime

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, os.path.abspath("_ext"))

import cellrank

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
needs_sphinx = "5.0"

# -- Project information -----------------------------------------------------

project = "CellRank"
author = cellrank.__author__
copyright = f"{datetime.now():%Y}, {author}"
release = cellrank.__version__
version = cellrank.__version__
github_repo = "cellrank"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx_gallery.load_style",
    "sphinx_design",
    "nbsphinx",
    "sphinx_copybutton",
    "typed_returns",
    "edit_on_github",
    "sphinxcontrib.bibtex",
]

intersphinx_mapping = {
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
    "squidpy": ("https://squidpy.readthedocs.io/en/latest/", None),
    "scvelo": ("https://scvelo.readthedocs.io/en/latest/", None),
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference/", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "statsmodels": ("https://www.statsmodels.org/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "joblib": ("https://joblib.readthedocs.io/en/latest/", None),
    "sklearn": ("https://scikit-learn.org/stable/", None),
    "seaborn": ("https://seaborn.pydata.org/", None),
    "pygam": ("https://pygam.readthedocs.io/en/latest/", None),
    "jax": ("https://jax.readthedocs.io/en/latest/", None),
    "pygpcca": ("https://pygpcca.readthedocs.io/en/latest/", None),
    "ot": ("https://pythonot.github.io/", None),
    "moscot": ("https://moscot.readthedocs.io/en/latest/", None),
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
source_suffix = [".rst", ".ipynb"]
master_doc = "index"

# syntax highlight
pygments_style = "default"
pygments_dark_style = "native"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["release/changelog/*", "_build", "**.ipynb_checkpoints"]

# bibliography
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"
bibtex_default_style = "alpha"

# linkcheck
user_agent = "Mozilla/5.0 (X11; Linux x86_64; rv:25.0) Gecko/20100101 Firefox/25.0"
# Twitter (used for handles in contributors.rst) doesn't like the above user-agent
linkcheck_ignore = [r"https://twitter\.com/.*", r"https://mobile\.twitter\.com/.*"]

# autodoc
autosummary_generate = True
autodoc_member_order = "bysource"
typehints_fully_qualified = False
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False

# theme
html_theme = "furo"
html_static_path = ["_static"]
html_logo = "_static/img/logo.png"
html_css_files = [
    "css/override.css",
    "css/sphinx_gallery.css",
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.1.1/css/all.min.css",
]
html_show_sphinx = False
html_show_sourcelink = True
html_theme_options = {
    "sidebar_hide_name": True,
    "light_css_variables": {
        "color-brand-primary": "#003262",
        "color-brand-content": "#003262",
        "admonition-font-size": "var(--font-size-normal)",
        "admonition-title-font-size": "var(--font-size-normal)",
        "code-font-size": "var(--font-size--small)",
    },
}

# spelling
spelling_lang = "en_US"
spelling_warning = True
spelling_word_list_filename = ["spelling/general.txt", "spelling/autosummary.txt"]
spelling_add_pypi_package_names = True
spelling_show_suggestions = False
spelling_exclude_patterns = ["references.rst"]
# see: https://pyenchant.github.io/pyenchant/api/enchant.tokenize.html
spelling_filters = [
    "enchant.tokenize.URLFilter",
    "enchant.tokenize.EmailFilter",
    "enchant.tokenize.MentionFilter",
    "docs.source.utils.ModnameFilter",
]

# nbsphinx
nbsphinx_execute = "never"
# TODO(michalk8): check this
nbsphinx_prolog = r"""
.. raw:: html
{{% set docname = env.doc2path(env.docname, base=None).split("/")[-1] %}}
.. raw:: html
    <style>
        p {{
            margin-bottom: 0.5rem;
        }}
        /* Main index page overview cards */
        /* https://github.com/spatialaudio/nbsphinx/pull/635/files */
        .jp-RenderedHTMLCommon table,
        div.rendered_html table {{
        border: none;
        border-collapse: collapse;
        border-spacing: 0;
        font-size: 12px;
        table-layout: fixed;
        color: inherit;
        }}
        body:not([data-theme=light]) .jp-RenderedHTMLCommon tbody tr:nth-child(odd),
        body:not([data-theme=light]) div.rendered_html tbody tr:nth-child(odd) {{
        background: rgba(255, 255, 255, .1);
        }}
    </style>
.. raw:: html
    <div class="admonition note">
        <p class="admonition-title">Note</p>
        <p>
        This page was generated from
        <a class="reference external" href="https://github.com/theislab/cellrank_notebooks/tree/{version}/">{docname}</a>.
        Interactive online version:
        <span style="white-space: nowrap;"><a href="https://colab.research.google.com/github/theislab/cellrank_notebooks/blob/{version}/{docname}"><img alt="Colab badge" src="https://colab.research.google.com/assets/colab-badge.svg" style="vertical-align:text-bottom"></a>.</span>
        Some tutorial content may look better in light mode.
        </p>
    </div>
""".format(
    version=version, docname="{{ docname|e }}"
)
