from typing import Any, Dict, List, Tuple, Union, Callable, Iterable, Optional

import wrapt
from types import MappingProxyType
from functools import wraps, update_wrapper
from contextlib import contextmanager
from multiprocessing import cpu_count

from anndata import AnnData
from cellrank import logging as logg
from anndata.utils import make_index_unique
from cellrank.ul._docs import d

import numpy as np
from scipy.sparse import issparse, spmatrix


def _check_collection(
    adata: AnnData,
    needles: Iterable[str],
    attr_name: str,
    key_name: str = "Gene",
    use_raw: bool = False,
    raise_exc: bool = True,
) -> List[str]:
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

    haystack, res = getattr(adata, attr_name), []
    for needle in needles:
        if needle not in haystack:
            if raise_exc:
                raise KeyError(
                    f"{key_name} `{needle}` not found in `{adata_name}.{attr_name}`."
                )
        else:
            res.append(needle)

    return res


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
    Minimum and maximum values, respectively.
    """
    if perc is not None:
        data = np.clip(data, *np.percentile(data, sorted(perc)))

    return float(np.nanmin(data)), float(np.nanmax(data))


def _modify_neigh_key(key: Optional[str]) -> str:
    if key in (None, "connectivities", "distances"):
        return "neighbors"

    if key.endswith("_connectivities"):
        key = key[:-15]
    elif key.endswith("_distances"):
        key = key[:-10]

    return key


def _get_neighs(
    adata: AnnData, mode: str = "distances", key: Optional[str] = None
) -> Union[np.ndarray, spmatrix]:
    if key is None:
        res = _read_graph_data(adata, mode)  # legacy behavior
    else:
        try:
            res = _read_graph_data(adata, key)
            assert isinstance(res, (np.ndarray, spmatrix))
        except (KeyError, AssertionError):
            res = _read_graph_data(adata, f"{_modify_neigh_key(key)}_{mode}")

    if not isinstance(res, (np.ndarray, spmatrix)):
        raise TypeError(
            f"Expected to find `numpy.ndarray` or `scipy.sparse.spmatrix`, found `{type(res)}`."
        )

    return res


def _has_neighs(adata: AnnData, key: Optional[str] = None) -> bool:
    return _modify_neigh_key(key) in adata.uns.keys()


def _get_neighs_params(adata: AnnData, key: str = "neighbors") -> Dict[str, Any]:
    return adata.uns.get(key, {}).get("params", {})


def _read_graph_data(adata: AnnData, key: str) -> Union[np.ndarray, spmatrix]:
    """
    Read graph data from :mod:`anndata`.

    Parameters
    ----------
    adata
        Annotated data object.
    key
        Key in :attr:`anndata.AnnData.obsp`.

    Returns
    -------
    The graph data.
    """
    if key in adata.obsp:
        return adata.obsp[key]

    raise KeyError(f"Unable to find data in `adata.obsp[{key!r}]`.")


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


@contextmanager
@d.dedent
def _gene_symbols_ctx(
    adata: AnnData,
    *,
    key: Optional[str] = None,
    use_raw: bool = False,
    make_unique: bool = False,
) -> AnnData:
    """
    Set gene names from a column in :attr:`anndata.AnnData.var`.

    Parameters
    ----------
    %(adata)s
    key
        Key in :attr:`anndata.AnnData.var` where the gene symbols are stored. If `None`, this operation is a no-op.
    use_raw
        Whether to change the gene names in :attr:`anndata.AnnData.raw`.
    make_unique
        Whether to make the newly assigned gene names unique.
    Yields
    ------
    The same ``adata`` with modified :attr:`anndata.AnnData.var_names`, depending on ``use_raw``.
    """

    def key_present() -> bool:
        if use_raw:
            if adata.raw is None:
                raise AttributeError(
                    "No `.raw` attribute found. Try specifying `use_raw=False`."
                )
            return key in adata.raw.var
        return key in adata.var

    if key is None:
        yield adata
    elif not key_present():
        raise KeyError(
            f"Unable to find gene symbols in `adata.{'raw.' if use_raw else ''}var[{key!r}]`."
        )
    else:
        adata_orig = adata
        if use_raw:
            adata = adata.raw

        var_names = adata.var_names.copy()
        try:
            # TODO(michalk8): doesn't update varm (niche)
            adata.var.index = (
                make_index_unique(adata.var[key]) if make_unique else adata.var[key]
            )
            yield adata_orig
        finally:
            # in principle we assume the callee doesn't change the index
            # otherwise, would need to check whether it has been changed and add an option to determine what to do
            adata.var.index = var_names


@wrapt.decorator
def _genesymbols(
    wrapped: Callable[..., Any], instance: Any, args: Any, kwargs: Any
) -> Any:
    with _gene_symbols_ctx(
        args[0],
        key=kwargs.pop("gene_symbols", None),
        make_unique=True,
        use_raw=kwargs.get("use_raw", False),
    ):
        return wrapped(*args, **kwargs)
