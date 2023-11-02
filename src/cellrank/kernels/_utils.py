from typing import Any, Callable

import wrapt

import numba as nb
import numpy as np
import pandas as pd
from numba import prange
from pandas.api.types import infer_dtype

from anndata import AnnData

from cellrank import logging as logg

jit_kwargs = {"nogil": True, "cache": True, "fastmath": True}


@nb.njit(parallel=False, **jit_kwargs)
def _np_apply_along_axis(func1d, axis: int, arr: np.ndarray) -> np.ndarray:
    """Apply a reduction function over a given axis.

    Parameters
    ----------
    func1d
        Reduction function that operates only on 1 dimension.
    axis
        Axis over which to apply the reduction.
    arr
        The array to be reduced.

    Returns
    -------
    The reduced array.
    """
    assert arr.ndim == 2
    assert axis in [0, 1]

    if axis == 0:
        result = np.empty(arr.shape[1])
        for i in range(len(result)):
            result[i] = func1d(arr[:, i])
        return result

    result = np.empty(arr.shape[0])
    for i in range(len(result)):
        result[i] = func1d(arr[i, :])

    return result


@nb.njit(**jit_kwargs)
def np_mean(array: np.ndarray, axis: int) -> np.ndarray:  # noqa
    return _np_apply_along_axis(np.mean, axis, array)


@nb.njit(**jit_kwargs)
def np_std(array: np.ndarray, axis: int) -> np.ndarray:  # noqa
    return _np_apply_along_axis(np.std, axis, array)


@nb.njit(**jit_kwargs)
def norm(array: np.ndarray, axis: int) -> np.ndarray:  # noqa
    return _np_apply_along_axis(np.linalg.norm, axis, array)


# this is faster than using flat array
@nb.njit(parallel=True)
def _random_normal(
    m: np.ndarray,
    v: np.ndarray,
    n_samples: int = 1,
) -> np.ndarray:
    """Sample number from normal distribution.

    Parameters
    ----------
    m
        Mean vector.
    v
        Variance vector.
    n_samples
        Number of samples to be generated.

    Returns
    -------
    Array of shape ``(n_samples x m.shape[0])``.
    """
    assert m.ndim == 1, "Means are not 1-dimensional."
    assert m.shape == v.shape, "Means and variances have different shape."

    if n_samples == 1:
        return np.expand_dims(np.array([np.random.normal(m[i], v[i]) for i in prange(m.shape[0])]), 0)  # noqa: NPY002

    return np.array(
        [[np.random.normal(m[i], v[i]) for _ in prange(n_samples)] for i in prange(m.shape[0])]  # noqa: NPY002
    ).T


@nb.njit(**jit_kwargs)
def _calculate_starts(indptr: np.ndarray, ixs: np.ndarray) -> np.ndarray:
    """Get the position where to put the data.

    Parameters
    ----------
    indptr
        Pointer of indices from :class:`~scipy.sparse.csr_matrix`.
    ixs
        Row indices for which to calculate the starts.

    Returns
    -------
    The starting positions.
    """
    starts = np.cumsum(indptr[ixs + 1] - indptr[ixs])
    return np.hstack((np.array([0], dtype=starts.dtype), starts))


def _get_basis(adata: AnnData, basis: str) -> np.ndarray:
    try:
        return adata.obsm[f"X_{basis}"]
    except KeyError:
        try:
            return adata.obsm[basis]  # e.g. 'spatial'
        except KeyError:
            raise KeyError(f"Unable to find a basis in `adata.obsm['X_{basis}']` or `adata.obsm[{basis!r}]`.") from None


def _ensure_numeric_ordered(adata: AnnData, key: str) -> pd.Series:
    if key not in adata.obs.keys():
        raise KeyError(f"Unable to find data in `adata.obs[{key!r}]`.")

    exp_time = adata.obs[key].copy()
    if not np.issubdtype(np.asarray(exp_time).dtype, np.number):
        try:
            exp_time = np.asarray(exp_time).astype(float)
        except Exception as e:  # noqa: BLE001/Cannot interpret
            raise TypeError(
                f"Unable to convert `adata.obs[{key!r}]` of type `{infer_dtype(adata.obs[key])}` to `float`."
            ) from e

    if not isinstance(exp_time.dtype, pd.CategoricalDtype):
        logg.debug(f"Converting `adata.obs[{key!r}]` to `categorical`")
        exp_time = np.asarray(exp_time)
        categories = sorted(set(exp_time[~np.isnan(exp_time)]))
        if len(categories) > 100:  # arbitrary threshold
            raise ValueError(
                f"Unable to convert `adata.obs[{key!r}]` to `categorical` since it "
                f"would create `{len(categories)}` categories."
            )
        exp_time = pd.Series(
            pd.Categorical(
                exp_time,
                categories=categories,
                ordered=True,
            )
        )

    if not exp_time.cat.ordered:
        logg.warning("Categories are not ordered. Using ascending order")
        exp_time.cat = exp_time.cat.as_ordered()

    exp_time = pd.Series(pd.Categorical(exp_time, ordered=True), index=adata.obs_names)
    if exp_time.isnull().any():
        raise ValueError("Series contains NaN value(s).")

    n_cats = len(exp_time.cat.categories)
    if n_cats < 2:
        raise ValueError(f"Expected to find at least `2` categories, found `{n_cats}`.")

    return exp_time


@wrapt.decorator
def require_tmat(
    wrapped: Callable[..., Any],
    instance: "KernelExpression",  # noqa: F821
    args: Any,
    kwargs: Any,
) -> Any:
    """Require that the transition matrix is computed before calling the wrapped function."""
    # this can trigger combinations, but not individual kernels
    if instance.transition_matrix is None:
        raise RuntimeError("Compute transition matrix first as `.compute_transition_matrix()`.")
    return wrapped(*args, **kwargs)
