# -*- coding: utf-8 -*-
"""General utility functions module."""
import pickle
from types import MappingProxyType
from typing import Any, Dict, Tuple, Union, TypeVar, Iterable, Optional
from pathlib import Path
from functools import wraps, update_wrapper
from multiprocessing import cpu_count

import numpy as np
from scipy.sparse import issparse, spmatrix

from cellrank import logging as logg
from cellrank.ul._docs import d

AnnData = TypeVar("AnnData")


class Pickleable:
    """Class which allows serialization and deserialization using :mod:pickle."""

    @d.get_full_descriptionf("pickleable")
    @d.get_sectionsf("pickleable", sections=["Parameters", "Returns"])
    def write(self, fname: Union[str, Path], ext: Optional[str] = "pickle") -> None:
        """
        Serialize self to a file.

        Parameters
        ----------
        fname
            Filename where to save the object.
        ext
            Filename extension to use. If `None`, don't append any extension.

        Returns
        -------
        None
            Nothing, just writes itself to a file using :mod:`pickle`.
        """

        fname = str(fname)
        if ext is not None:
            if not ext.startswith("."):
                ext = "." + ext
            if not fname.endswith(ext):
                fname += ext

        logg.debug(f"Writing to `{fname}`")

        with open(fname, "wb") as fout:
            pickle.dump(self, fout)

    @staticmethod
    def read(fname: Union[str, Path]) -> Any:
        """
        Deserialize self from a file.

        Parameters
        ----------
        fname
            Filename from which to read the object.

        Returns
        -------
        :class:`typing.Any`
            The deserialized object.
        """

        with open(fname, "rb") as fin:
            return pickle.load(fin)


def _check_collection(
    adata: AnnData,
    needles: Iterable[str],
    attr_name: str,
    key_name: str = "Gene",
    use_raw: bool = False,
) -> None:
    """
    Check if given collection contains all the keys.

    Parameters
    ----------
    adata: :class:`anndata.AnnData`
        Annotated data object.
    needles
        Keys to check.
    attr_name
        Attribute of ``adata`` where the needles are stored.
    key_name
        Pretty name of the key which will be displayed when error is found.
    use_raw
        Whether to access ``adata.raw`` or just ``adata``.

    Returns
    -------
    None
        Nothing, but raises and :class:`KeyError` if one of the needles is not found.
    """
    adata_name = "adata"

    if use_raw and adata.raw is None:
        logg.warning(
            "Argument `use_raw` was set to `True`, but no `raw` attribute is found. Ignoring"
        )
        use_raw = False
    if use_raw:
        adata_name = "adata.raw"
        adata = adata.raw

    haystack = getattr(adata, attr_name)
    for needle in needles:
        if needle not in haystack:
            raise KeyError(
                f"{key_name} `{needle}` not found in `{adata_name}.{attr_name}`."
            )


def _get_n_cores(n_cores: Optional[int], n_jobs: Optional[int]) -> int:
    """
    Make number of cores a positive integer.

    Parameters
    ----------
    n_cores
        Number of cores to use.
    n_jobs.
        Number of jobs. This is just used to determine if the collection is a singleton.
        If `1`, always returns `1`.

    Returns
    -------
    int
        Positive integer corresponding to how many cores to use.
    """
    if n_cores == 0:
        raise ValueError("Number of cores cannot be `0`.")
    if n_jobs == 1 or n_cores is None:
        return 1
    if n_cores < 0:
        return cpu_count() + 1 + n_cores

    return n_cores


def _minmax(
    data: np.ndarray, perc: Optional[Tuple[float, float]] = None
) -> Tuple[float, float]:
    """
    Return minimum and maximum value of the data.

    Parameters
    ----------
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


def _has_neighs(adata: AnnData) -> bool:
    return "neighbors" in adata.uns.keys()


def _get_neighs(adata: AnnData, mode: str = "distances") -> Union[np.ndarray, spmatrix]:
    try:
        return _read_graph_data(adata, mode)
    except KeyError:
        return _read_graph_data(adata, "neighors")[mode]


def _get_neighs_params(adata: AnnData) -> Dict[str, Any]:
    return adata.uns["neighbors"]["params"]


def _read_graph_data(adata: AnnData, key: str) -> Union[np.ndarray, spmatrix]:
    """
    Read graph data from :mod:`anndata`.

    :module`AnnData` >=0.7 stores `(n_obs x n_obs)` matrices in `.obsp` rather than `.uns`.
    This is for backward compatibility.

    Parameters
    ----------
    adata
        Annotated data object.
    key
        Key from either ``adata.uns`` or ``adata.obsp``.

    Returns
    -------
    :class:`numpy.ndarray` or :class:`scipy.sparse.spmatrix`
        The graph data.
    """

    if hasattr(adata, "obsp") and adata.obsp is not None and key in adata.obsp.keys():
        logg.debug(f"Reading key `{key!r}` from `adata.obsp`")
        return adata.obsp[key]

    if key in adata.uns.keys():
        logg.debug(f"Reading key `{key!r}` from `adata.uns`")
        return adata.uns[key]

    raise KeyError(f"Unable to find key `{key!r}` in `adata.obsp` or `adata.uns`.")


def _write_graph_data(
    adata: AnnData,
    data: Union[np.ndarray, spmatrix],
    key: str,
):
    """
    Write graph data to :mod:`AnnData`.

    :module`anndata` >=0.7 stores `(n_obs x n_obs)` matrices in `.obsp` rather than `.uns`.
    This is for backward compatibility.

    Parameters
    ----------
    adata
        Annotated data object.
    data
        The graph data we want to write.
    key
        Key from either ``adata.uns`` or `adata.obsp``.

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

    logg.debug(f"Writing graph data to `adata.{write_to}[{key!r}]`")


def valuedispatch(func):
    """Dispatch a function based on the first value."""
    registry = {}

    def dispatch(value):
        return registry.get(value, func)

    def register(value, func=None):
        if func is None:
            return lambda f: register(value, f)
        registry[value] = func
        return func

    @wraps(func)
    def wrapper(*args, **kwargs):
        return dispatch(args[0])(*args[1:], **kwargs)

    wrapper.register = register
    wrapper.dispatch = dispatch
    wrapper.registry = MappingProxyType(registry)

    update_wrapper(wrapper, func)

    return wrapper


def _densify_squeeze(x: Union[spmatrix, np.ndarray], dtype=np.float32) -> np.ndarray:
    if issparse(x):
        x = x.toarray()
    # use np.array instead of asarray to create a copy
    x = np.array(x, dtype=dtype)
    if x.ndim == 2 and x.shape[1] == 1:
        x = np.squeeze(x, axis=1)

    return x
