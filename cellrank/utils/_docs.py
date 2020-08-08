# -*- coding: utf-8 -*-
"""Module for documentation helper function."""

from textwrap import dedent

from docrep import DocstringProcessor

_adata = """\
adata : :class:`~anndata.AnnData`
    Annotated data object."""
_plotting = """\
figsize
    Size of the figure.
dpi
    Dots per inch.
save
    Filename where to save the plot."""
_parallel = """\
show_progress_bar
    Whether to show a progress bar. Disabling it may improve performance.
n_jobs
    Number of parallel jobs. If `-1`, use all available cores. If `None` or `1`, the execution is sequential.
backend
    Which backend to use for multiprocessing. See :class:`joblib.Parallel` for valid options."""
_model = """\
model
    BaseModel to fit.

    - If a :class:`dict`, gene and lineage specific models can be specified. Use `'*'` to indicate
      all genes or lineages, for example `{'Map2': {'*': ...}, 'Dcx': {'Alpha': ..., '*': ...}}`."""
_just_plots = """\
None
    Nothing, just plots the figure. Optionally saves it based on :paramref:`save`."""
_backward = """\
backward
    Direction of the process."""
_eigen = """\
which
    Eigenvalues are in general complex. `'LR'` - largest real part, `'LM'` - largest magnitude.
alpha
    Used to compute the `eigengap`. :paramref:`alpha` is the weight given
    to the deviation of an eigenvalue from one."""
_n_cells = """\
n_cells
    Number of most likely cells from each main to select."""
_fit = """\
n_lineages
    Number of lineages. If `None`, it will be determined automatically.
cluster_key
    Match computed states against pre-computed clusters to annotate the states.
    For this, provide a key from :paramref:`adata` `.obs` where cluster labels have been computed.
keys
    Determines which %(root_or_final) states to use by passing their names.
    Further, %(root_or_final)s states can be combined. If e.g. the %(final)s states are
    ['Neuronal_1', 'Neuronal_1', 'Astrocytes', 'OPC'], then passing keys=['Neuronal_1, Neuronal_2', 'OPC']
    means that the two neuronal %(final)s states are treated as one and the 'Astrocyte' state is excluded.
"""
_density_correction = (
    "Optionally, we apply a density correction as described in [Coifman05]_, "
    "where we use the implementation of [Haghverdi16]_."
)
_time_range = """
time_range
   If a :class:`tuple`, it specifies the minimum and maximum pseudotime. Both values can be `None`, in which case
   the minimum is the minimum pseudotime and maximum is automatically determined. If :class:`float`,
   it specified the maximum pseudotime.
"""

_copy = """Return a copy of self."""
_root = "root"
_final = "final"


def inject_docs(**kwargs):
    r"""Docstrings should start with "\" in the first line for proper formatting."""

    def decorator(obj):
        obj.__doc__ = dedent(obj.__doc__).format(**kwargs)
        return obj

    return decorator


d = DocstringProcessor(
    plotting=_plotting,
    parallel=_parallel,
    model=_model,
    adata=_adata,
    just_plots=_just_plots,
    backward=_backward,
    root=_root,
    final=_final,
    eigen=_eigen,
    root_or_final=f"{_root} or {_final}",
    n_cells=_n_cells,
    fit=_fit,
    copy=_copy,
    density_correction=_density_correction,
)
