# -*- coding: utf-8 -*-
from multiprocessing import cpu_count
from typing import Iterable, Hashable, Dict, Optional, Tuple, List, Union, Any

from scipy.sparse import spmatrix
import anndata
import numpy as np


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


def has_neighs(adata: anndata.AnnData) -> bool:
    return "neighbors" in adata.uns.keys()


def get_neighs(
    adata: anndata.AnnData, mode: str = "distances"
) -> Union[np.ndarray, spmatrix]:
    return (
        adata.obsp[mode]
        if hasattr(adata, "obsp") and mode in adata.obsp.keys()
        else adata.uns["neighbors"][mode]
    )


def get_neighs_params(adata: anndata.AnnData) -> Dict[str, Any]:
    return adata.uns["neighbors"]["params"]
