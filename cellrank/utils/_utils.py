# -*- coding: utf-8 -*-
"""General utility functions module."""

from typing import Any, Dict, List, Tuple, Union, Hashable, Iterable, Optional
from multiprocessing import cpu_count

import numpy as np
from scipy.sparse import spmatrix

import anndata
from scanpy import logging as logg
from anndata import AnnData


def check_collection(
    adata: anndata.AnnData,
    needles: Iterable[str],
    attr_name: str,
    key_name: str = "Gene",
) -> None:
    """
    Check if given collection contains all the keys.

    Params
    ------
    adata: :class:`anndata.AnnData`
        Annotated data object.
    needles
        Keys to check.
    attr_name
        Attribute of :paramref:`adata` where the needles are stored.

    Returns
    -------
    None
        Nothing, but raises and :class:`KeyError` if one of needles is not found.
    """

    haystack = getattr(adata, attr_name)
    for needle in needles:
        if needle not in haystack:
            raise KeyError(f"{key_name} `{needle}` not found in `adata.{attr_name}`.")


def _get_n_cores(n_cores: Optional[int], n_genes: int) -> int:
    """
    Make number of cores a positive integer.

    Params
    ------
    n_cores
        Number of cores to use.
    n_genes.
        Number of genes.

    Returns
    -------
    int
        Positive integer corresponding to how many cores to use.
    """

    if n_cores == 0:
        raise ValueError("Number of cores cannot be `0`.")
    if n_genes == 1 or n_cores is None:
        return 1
    if n_cores < 0:
        return cpu_count() + 1 + n_cores

    return n_cores


def _minmax(
    data: np.ndarray, perc: Optional[Tuple[float, float]] = None
) -> Tuple[float, float]:
    """
    Return minimum and maximum value of the data.

    Params
    ------
    data
        Values for which to return the minimum and maximum.
    perc
        If not `None`, clip the values by the percentiles before getting the result.

    Returns
    -------
    :class:`tuple`
        Minimum and maximum values, respectively.
    """

    if perc is not None:
        data = np.clip(data, *np.percentile(data, sorted(perc)))

    return float(np.nanmin(data)), float(np.nanmax(data))


def _make_unique(collection: Iterable[Hashable]) -> List[Hashable]:
    """
    Make a collection unique while maintaining the order.

    Params
    ------
    collection
        Values to make unique.

    Returns
    -------
    :class:`list`
        The same collection with unique values.
    """

    seen, res = set(), []
    for item in collection:
        if item not in seen:
            seen.add(item)
            res.append(item)

    return res


def _has_neighs(adata: anndata.AnnData) -> bool:
    return "neighbors" in adata.uns.keys()


def _get_neighs(
    adata: anndata.AnnData, mode: str = "distances"
) -> Union[np.ndarray, spmatrix]:
    try:
        return _read_graph_data(adata, mode)
    except KeyError:
        return _read_graph_data(adata, "neighors")[mode]


def _get_neighs_params(adata: anndata.AnnData) -> Dict[str, Any]:
    return adata.uns["neighbors"]["params"]


def _read_graph_data(adata: AnnData, key: str) -> Union[np.ndarray, spmatrix]:
    """
    Read graph data from :module:`anndata`.

    :module`AnnData` `>=0.7` stores (n_obs x n_obs) matrices in `.obsp` rather than `.uns`.
    This is for backward compatibility.

    Params
    ------
    adata
        Annotated data object.
    key
        Key from either `.uns` or `.obsp`.

    Returns
    -------
    :class:`numpy.ndarray` or :class:`scipy.sparse.spmatrix`
        The graph data.
    """

    if hasattr(adata, "obsp") and adata.obsp is not None and key in adata.obsp.keys():
        logg.debug(f"Read key `{key!r}` from `adata.obsp`")
        return adata.obsp[key]

    if key in adata.uns.keys():
        logg.debug(f"Read key `{key!r}` from `adata.uns`")
        return adata.uns[key]

    raise KeyError(f"Unable to find key `{key!r}` in `adata.obsp` or `adata.uns`.")


def _write_graph_data(
    adata: AnnData, data: Union[np.ndarray, spmatrix], key: str,
):
    """
    Write graph data to :module:`AnnData`.

    :module`anndata` `>=0.7` stores (n_obs x n_obs) matrices in `.obsp` rather than `.uns`.
    This is for backward compatibility.

    Params
    ------
    adata
        Annotated data object.
    data
        The graph data we want to write.
    key
        Key from either `.uns` or `.obsp`.

    Returns
    --------
    None
        Nothing, just writes the data.
    """

    try:
        adata.obsp[key] = data
        write_to = "obsp"

        if data.shape[0] != data.shape[1]:
            logg.warning(
                f"`adata.obsp` attribute should only contain square matrices, found shape `{data.shape}`"
            )

    except AttributeError:
        adata.uns[key] = data
        write_to = "uns"

    logg.debug(f"Write graph data {key!r} to `adata.{write_to}`")
