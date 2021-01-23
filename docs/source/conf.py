import os
import sys
import logging
from pathlib import Path
from datetime import datetime
from collections import ChainMap, defaultdict
from urllib.parse import urljoin
from urllib.request import urlretrieve

from sphinx_gallery.sorting import ExplicitOrder, _SortKey

HERE = Path(__file__).parent
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
author = cellrank.__author__
copyright = f"{datetime.now():%Y}, {author}."
release = "master"
version = f"master ({cellrank.__version__})"


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
    "sphinx_paramlinks",
    "sphinx.ext.autosummary",
    "sphinx_gallery.gen_gallery",
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
    networkx=("https://networkx.org/documentation/stable/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    statsmodels=("https://www.statsmodels.org/stable/", None),
    matplotlib=("https://matplotlib.org/", None),
    joblib=("https://joblib.readthedocs.io/en/latest/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
    seaborn=("https://seaborn.pydata.org/", None),
    pygam=("https://pygam.readthedocs.io/en/latest/", None),
    jax=("https://jax.readthedocs.io/en/latest/", None),
    pygpcca=("https://pygpcca.readthedocs.io/en/latest/", None),
)

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates", "_build"]
source_suffix = [".rst", ".ipynb"]
master_doc = "index"
pygments_style = "sphinx"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "**.md5",
    "**.py",
    "**.ipynb_checkpoints",
]  # ignore anything that isn't .rst
suppress_warnings = ["ref.citation"]

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

# -- sphinx gallery


def reset_scvelo(gallery_conf, fname):
    import scvelo as scv

    scv.set_figure_params(
        style="scvelo", color_map="viridis", format="png", ipython_format="png"
    )


def reset_matplotlib(gallery_conf, fname):
    import matplotlib as mpl

    mpl.use("agg")

    import matplotlib.pyplot as plt

    plt.rcdefaults()
    mpl.rcParams["savefig.bbox"] = "tight"
    mpl.rcParams["savefig.transparent"] = True


example_dir = HERE.parent.parent / "examples"
rel_example_dir = Path("..") / ".." / "examples"


class ExplicitSubsectionOrder(_SortKey):

    _order = ChainMap(
        {
            example_dir / "estimators" / "compute_eigendecomposition.py": 0,
            example_dir / "estimators" / "compute_schur_vectors.py": 10,
            example_dir / "estimators" / "compute_schur_matrix.py": 20,
            example_dir / "estimators" / "compute_macrostates.py": 30,
            example_dir / "estimators" / "compute_coarse_T.py": 40,
            example_dir / "estimators" / "compute_terminal_states_gpcca.py": 50,
            example_dir / "estimators" / "compute_abs_probs.py": 60,
            example_dir / "estimators" / "compute_lineage_drivers.py": 70,
            example_dir / "estimators" / "compute_fit.py": 80,
        },
        {
            example_dir / "plotting" / "plot_initial_states.py": 0,
            example_dir / "plotting" / "plot_terminal_states.py": 10,
            example_dir / "plotting" / "plot_lineages.py": 20,
            example_dir / "plotting" / "plot_lineage_drivers.py": 30,
            example_dir / "plotting" / "plot_directed_paga.py": 40,
            example_dir / "plotting" / "plot_gene_trends.py": 50,
            example_dir / "plotting" / "plot_heatmap.py": 60,
            example_dir / "plotting" / "plot_cluster_lineage.py": 70,
            example_dir / "plotting" / "plot_cluster_fates.py": 80,
            example_dir / "plotting" / "plot_graph.py": 90,
        },
        defaultdict(
            lambda: 1000,
            {
                example_dir / "other" / "plot_model.py": 0,
                example_dir / "other" / "compute_lineage_tricks.py": 10,
                example_dir / "other" / "compute_kernel_tricks.py": 20,
            },
        ),
    )

    def __call__(self, filename: str) -> int:
        src_file = os.path.normpath(os.path.join(self.src_dir, filename))
        return self._order[Path(src_file)]

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {repr(dict(self._order))}>"


sphinx_gallery_conf = {
    "image_scrapers": "matplotlib",
    "reset_modules": (
        "seaborn",
        reset_matplotlib,
        reset_scvelo,
    ),  # reset scvelo as last
    "filename_pattern": f"{os.path.sep}(plot_|compute_)",
    "examples_dirs": example_dir,
    "gallery_dirs": "auto_examples",  # path to where to save gallery generated output
    "abort_on_example_error": True,
    "show_memory": True,
    "within_subsection_order": ExplicitSubsectionOrder,
    "subsection_order": ExplicitOrder(
        [
            rel_example_dir / "estimators",  # really must be relative
            rel_example_dir / "plotting",
            rel_example_dir / "other",
        ]
    ),
    "reference_url": {
        "sphinx_gallery": None,
    },
    "line_numbers": False,
    "compress_images": ("images", "thumbnails", "-o3"),
    "inspect_global_variables": False,
    "backreferences_dir": "gen_modules/backreferences",
    "doc_module": "cellrank",
    "download_all_examples": False,
    "pypandoc": True,  # convert rST to md when downloading notebooks
}

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
