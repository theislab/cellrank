# -*- coding: utf-8 -*-
"""Module for documentation helper function."""

from textwrap import dedent

from docrep import DocstringProcessor

_adata = """\
adata : :class:`anndata.AnnData`
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
    Whether to show a progress bar when parallelizing.
n_jobs
    Number of parallel jobs. If `-1`, use all available cores. If `None` or `1`, the execution is sequential.
backend
    Which backend to use for multiprocessing. See :class:`joblib.Parallel` for valid options."""
_model = """\
model
    Model to fit.

    - If a :class:`dict`, gene and lineage specific models can be specified. Use `'*'` to indicate
    all genes or lineages, for example `{'Map2': {'*': ...}, 'Dcx': {'Alpha': ..., '*': ...}}`."""
_just_plots = """\
None
    Nothing, just plots the figure. Optionally saves it based on :paramref:`save`."""
_final = """\
final
    Whether to consider cells going to final states or vice versa."""


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
    final=_final,
)
