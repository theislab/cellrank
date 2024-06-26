# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import sys
from datetime import datetime
from pathlib import Path

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
import cellrank

sys.path.insert(0, str(Path(__file__).parent / "_ext"))

# -- Project information -----------------------------------------------------

project = cellrank.__name__
author = cellrank.__author__
version = cellrank.__version__
copyright = f"{datetime.now():%Y}, Theislab"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinxcontrib.bibtex",
    "sphinx_copybutton",
    "sphinx_autodoc_typehints",
    "myst_nb",
    "sphinx_tippy",
    "sphinx_design",
    "typed_returns",
]
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "jax": ("https://jax.readthedocs.io/en/latest/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "sklearn": ("https://scikit-learn.org/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "statsmodels": ("https://www.statsmodels.org/stable/", None),
    "pygam": ("https://pygam.readthedocs.io/en/latest/", None),
    "pygpcca": ("https://pygpcca.readthedocs.io/en/latest/", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
    "joblib": ("https://joblib.readthedocs.io/en/latest/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "seaborn": ("https://seaborn.pydata.org/", None),
    "anndata": ("https://anndata.readthedocs.io/en/latest/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/latest/", None),
    "scvelo": ("https://scvelo.readthedocs.io/en/latest/", None),
    "squidpy": ("https://squidpy.readthedocs.io/en/latest/", None),
    "moscot": ("https://moscot.readthedocs.io/en/latest/", None),
    "ot": ("https://pythonot.github.io/", None),
}
master_doc = "index"
pygments_style = "tango"
pygments_dark_style = "monokai"

nitpicky = True

# bibliography
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"
bibtex_default_style = "alpha"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
}

# myst
nb_execution_mode = "off"
myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "amsmath",
]
myst_heading_anchors = 2

# hover
tippy_anchor_parent_selector = "div.content"
tippy_enable_mathjax = True
# no need because of sphinxcontrib-bibtex
tippy_enable_doitips = False

# autodoc + napoleon
autosummary_generate = True
autodoc_member_order = "alphabetical"
autodoc_typehints = "description"
autodoc_mock_imports = ["moscot"]
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# spelling
spelling_lang = "en_US"
spelling_warning = True
spelling_word_list_filename = "spelling_wordlist.txt"
spelling_add_pypi_package_names = True
spelling_exclude_patterns = ["references.rst"]
# see: https://pyenchant.github.io/pyenchant/api/enchant.tokenize.html
spelling_filters = [
    "enchant.tokenize.URLFilter",
    "enchant.tokenize.EmailFilter",
    "enchant.tokenize.MentionFilter",
]

linkcheck_ignore = [
    # 403 Client Error: Forbidden for url
    r"https://doi.org/10.1063/1.5064530",
    r"https://doi.org/10.1021/acs.jctc.8b00079",
    r"https://doi.org/10.1093/bioinformatics/btv325",
    r"https://doi.org/10.1242/dev.173849",
    r"https://doi.org/10.1242/dev.126011",
    r"https://doi.org/10.2337/db20-0599",
    r"https://doi.org/10.1210/en.2013-1663",
    r"https://doi.org/10.4161/isl.21984",
    r"https://doi.org/10.15252/msb.202110282",
    r"https://doi.org/10.1073/pnas.0500334102",
    r"https://doi.org/10.1126/science.aar3131",
    r"https://doi.org/10.1126/science.aax0249",
    r"https://doi.org/10.1126/science.aax3072",
    r"https://www.science.org/doi/full/10.1126/science.1247651",
]

exclude_patterns = [
    "_build",
    "notebooks/README.rst",
    "notebooks/notebooks/MK_2020-10-19_generate_docs_figure.ipynb",
    "release/changelog/*",
    "**.ipynb_checkpoints",
]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "furo"
html_static_path = ["_static"]
html_css_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.1.1/css/all.min.css",
    "css/override.css",
]

html_show_sphinx = False
html_show_sourcelink = False
html_theme_options = {
    "announcement": "If you're moving from CellRank v1 to v2, "
    "see our <a href='https://cellrank.readthedocs.io/en/latest/about/version2.html'>notes on important changes</a>.",
    "sidebar_hide_name": True,
    "light_logo": "img/light_mode_logo.png",
    "dark_logo": "img/dark_mode_logo.png",
    "light_css_variables": {
        "color-brand-primary": "#003262",
        "color-brand-content": "#003262",
        "admonition-font-size": "var(--font-size-normal)",
        "admonition-title-font-size": "var(--font-size-normal)",
    },
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/theislab/cellrank",
            "html": "",
            "class": "fab fa-github",
        },
    ],
}
