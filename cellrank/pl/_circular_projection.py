# -*- coding: utf-8 -*-
from types import MappingProxyType
from typing import Any, Tuple, Union, Mapping, Callable, Optional, Sequence
from pathlib import Path

import scvelo as scv
from anndata import AnnData

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances

import matplotlib.pyplot as plt

from cellrank import logging as logg
from cellrank.tl import Lineage
from cellrank.ul._docs import d
from cellrank.pl._utils import _held_karp
from cellrank.tl._utils import save_fig, _unique_order_preserving
from cellrank.ul._utils import _check_collection
from cellrank.tl._constants import ModeEnum, AbsProbKey


class LineageOrder(ModeEnum):  # noqa: D101
    DEFAULT = "default"
    OPTIMAL = "optimal"


class SpecialKey(ModeEnum):  # noqa: D101
    DEGREE = "priming_degree"
    DIRECTION = "priming_direction"


Metric_T = Union[str, Callable, np.ndarray, pd.DataFrame]


def _get_distances(data: Lineage, metric: Metric_T) -> np.ndarray:
    if isinstance(metric, str) or callable(metric):
        metric = pairwise_distances(data.X.T, metric=metric)
    elif isinstance(metric, np.ndarray):
        shape = (data.shape[1], data.shape[1])
        if metric.shape != shape:
            raise ValueError(
                f"Expected a `numpy.ndarray` of shape `{shape}`, found `{metric.shape}`."
            )
    elif isinstance(metric, pd.DataFrame):
        if np.any(metric.index != data.names):
            raise ValueError(
                f"Expected `pandas.DataFrame` to have the following index `{list(data.names)}`, "
                f"found `{list(metric.columns)}`."
            )
        if np.any(metric.columns != data.names):
            raise ValueError(
                f"Expected `pandas.DataFrame` to have the following columns `{list(data.names)}`, "
                f"found `{list(metric.columns)}`."
            )
    else:
        raise TypeError(
            f"Expected either metric defined by `str`, `callable` or a pairwise distance matrix of type"
            f" `numpy.ndarray` or `pandas.DataFrame`, found `{type(metric).__name__}`."
        )

    return np.asarray(metric, dtype=np.float64)


def _get_optimal_order(data: Lineage, metric: Metric_T) -> Tuple[float, np.ndarray]:
    """Solve the TSP using dynamic programming."""
    return _held_karp(_get_distances(data, metric))


def _priming_degree(probs: Union[np.ndarray, Lineage]) -> np.ndarray:
    if isinstance(probs, Lineage):
        probs = probs.X
    return -np.sum(probs * np.log(probs / np.mean(probs, axis=0)), axis=1)


# TODO: do we event want this?
def _priming_direction(probs: Union[np.ndarray, Lineage]) -> np.ndarray:
    if isinstance(probs, Lineage):
        probs = probs.X
    return np.argmax(probs / np.sum(probs, axis=0), axis=1)


@d.dedent
def circular_projection(
    adata: AnnData,
    keys: Union[str, Sequence[str]],
    backward: bool = False,
    delta: Optional[float] = 1.25,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    lineage_order: Optional[str] = None,
    metric: Union[str, Callable, np.ndarray, pd.DataFrame] = "correlation",
    ncols: int = 4,
    space: float = 0.25,
    use_raw: bool = False,
    text_kwargs: Mapping[str, Any] = MappingProxyType({}),
    key_added: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, Path]] = None,
    **kwargs,
):
    r"""
    Plot absorption probabilities on a circular embedding.

    Parameters
    ----------
    %(adata)s
    keys
        Keys in :attr:`anndata.AnnData.obs` or :attr:`anndata.AnnData.var_names`. Extra available keys:

            - `'priming_degree'`: the entropy relative to the "average" cell, as described in [Velten17]_. It follows \
            the formula :math:`\sum_{j} p_{i, j} \log \frac{p_{ij}}{\overline{p_j}}` where :math:`i` corresponds to a \
            cell and :math:`j` corresponds to a lineage and :math:`\overline{p_j}` is the average probability for
            lineage :math:`j`.
    %(backward)s
    delta
        Distance at which the lineage labels will be drawn. If `None`, don't draw the lineage labels.
    lineages
        Lineages to plot. If `None`, plot all lineages.
    lineage_order
        Can be one of the following:

            - `None`: it will determined automatically, based on the number of lineages.
            - `'optimal'`: order the lineages optimally by solving the Travelling salesman problem. \
            Recommended for <= `20` lineages.
            - `'default'`: use the order as specified in ``lineages``.
    metric
        Metric to use when contructing pairwise distance matrix when ``lineage_order = 'optimal'``. For available
        options, see :func:`sklearn.metrics.pairwise_distances`.
    ncols
        Number of columns when plotting multiple ``keys``.
    space
        Horizontal and vertical space between for :func:`matplotlib.pyplot.subplots_adjust`.
    use_raw
        Whether to access :attr:`anndata.AnnData.raw` when there are ``keys`` in :attr:`anndata.AnnData.var_names`.
    text_kwargs
        Keyword arguments for :func:`matplotlib.pyplot.text`.
    key_added
        Key in :attr:`anndata.AnnData.obsm` where to add the circular embedding. If `None`, it will be set to
        `'X_fate_simplex_{fwd,bwd}'`, based on ``backward``.
    %(plotting)s
    **kwargs
        Keyword arguments for :func:`scvelo.pl.scatter`.

    Returns
    -------
    %(just_plots)s
    """
    # TODO: logging, tests, polish docstrings, references
    # TODO: warn if using optimal and too many lineages (> 10)
    if key_added is None:
        key_added = "X_fate_simplex_" + ("bwd" if backward else "fwd")

    if isinstance(keys, str):
        keys = (keys,)

    keys = _unique_order_preserving(keys)
    keys_ = _check_collection(
        adata, keys, "obs", key_name="Observation", raise_exc=False
    ) + _check_collection(
        adata, keys, "var_names", key_name="Gene", raise_exc=False, use_raw=use_raw
    )
    haystack = {s.s for s in SpecialKey}
    keys = keys_ + [k for k in keys if k in haystack]

    if not len(keys):
        raise ValueError("No valid keys have been selected.")

    lineage_key = str(AbsProbKey.BACKWARD if backward else AbsProbKey.FORWARD)
    if lineage_key not in adata.obsm:
        raise KeyError(f"Lineages key `{lineage_key!r}` not found in `adata.obsm`.")

    probs = adata.obsm[lineage_key]

    if isinstance(lineages, str):
        lineages = (lineages,)
    elif lineages is None:
        lineages = probs.names

    probs = adata.obsm[lineage_key][lineages]
    n_lin = probs.shape[1]
    if n_lin <= 1:
        raise ValueError(f"Expected at least `2` lineages, found `{n_lin}`")

    if lineage_order is None:
        lineage_order = LineageOrder.OPTIMAL if n_lin <= 15 else LineageOrder.DEFAULT
    lineage_order = LineageOrder(lineage_order)

    if lineage_order == LineageOrder.OPTIMAL:
        logg.info(f"Solving TSP for `{n_lin}` states")
        _, order = _get_optimal_order(probs, metric=metric)
    else:
        order = np.arange(n_lin)

    probs = probs[:, order]

    angle_vec = np.linspace(0, 2 * np.pi, n_lin, endpoint=False)
    angle_vec_sin = np.sin(angle_vec)
    angle_vec_cos = np.cos(angle_vec)

    x = np.sum(probs.X * angle_vec_sin, axis=1)
    y = np.sum(probs.X * angle_vec_cos, axis=1)
    adata.obsm[key_added] = np.c_[x, y]

    nrows = int(np.ceil(len(keys) / ncols))
    fig, ax = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(ncols * 5, nrows * 5) if figsize is None else figsize,
        dpi=dpi,
    )

    fig.subplots_adjust(wspace=space, hspace=space)
    axes = np.ravel([ax])

    _i = 0
    for _i, (k, ax) in enumerate(zip(keys, axes)):

        try:
            k = SpecialKey(k)
            kwargs["title"] = k
            if k == SpecialKey.DEGREE:
                k = _priming_degree(probs)
            else:
                # TODO: for SK.Direction, how to supply to the categorical colors, should we just save to anndata?
                raise NotImplementedError(k)
        except ValueError:
            # not a special key
            pass

        scv.pl.scatter(
            adata,
            basis=key_added,
            color=k,
            show=False,
            ax=ax,
            use_raw=use_raw,
            **kwargs,
        )

        if delta is None:
            continue
        for j, (lin, c) in enumerate(zip(probs.names, probs.colors)):
            ax.text(
                delta * angle_vec_sin[j],
                delta * angle_vec_cos[j],
                s=lin,
                ha="center",
                c=c,
                **text_kwargs,
            )

    for j in range(_i + 1, len(axes)):
        axes[j].remove()

    if save is not None:
        save_fig(fig, save)
