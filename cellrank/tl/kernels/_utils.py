# -*- coding: utf-8 -*-
"""Utility functions for kernels, mostly VelocityKernel."""
from typing import Tuple, Callable, Optional
from inspect import signature

import numpy as np
from numba import njit, prange
from scipy.sparse import csr_matrix

jit_kwargs = {"nogil": True, "cache": True, "fastmath": True}


try:
    import jax.numpy as jnp
    from jax import jit, hessian

    _HAS_JAX = True

    @jit
    def _predict_transition_probabilities_jax(
        X: np.ndarray,
        W: np.ndarray,
        softmax_scale: float = 1,
    ):
        # pearson correlation
        W -= W.mean(axis=1)[:, None]
        X -= X.mean()

        W_norm = jnp.linalg.norm(W, axis=1)
        X_norm = jnp.linalg.norm(X)
        denom = X_norm * W_norm

        mask = jnp.isclose(denom, 0)
        denom = jnp.where(jnp.isclose(denom, 0), 1, denom)  # essential

        x = W.dot(X) / denom

        numerator = x * softmax_scale
        numerator = jnp.exp(numerator - jnp.nanmax(numerator))
        numerator = jnp.where(mask, 0, numerator)  # essential

        softmax = numerator / jnp.nansum(numerator)

        return softmax

    _predict_transition_probabilities_jax_H = hessian(
        _predict_transition_probabilities_jax
    )

except ImportError:
    _HAS_JAX = False
    _predict_transition_probabilities_jax = None
    _predict_transition_probabilities_jax_H = None


@njit(parallel=False, **jit_kwargs)
def _np_apply_along_axis(func1d, axis: int, arr: np.ndarray) -> np.ndarray:
    """
    Apply a reduction function over a given axis.

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
    :class:`numpy.ndarray`
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


@njit(**jit_kwargs)
def np_mean(array: np.ndarray, axis: int) -> np.ndarray:  # noqa
    return _np_apply_along_axis(np.mean, axis, array)


@njit(**jit_kwargs)
def np_max(array: np.ndarray, axis: int) -> np.ndarray:  # noqa
    return _np_apply_along_axis(np.max, axis, array)


@njit(**jit_kwargs)
def np_sum(array: np.ndarray, axis: int) -> np.ndarray:  # noqa
    return _np_apply_along_axis(np.sum, axis, array)


@njit(**jit_kwargs)
def norm(array: np.ndarray, axis: int) -> np.ndarray:  # noqa
    return _np_apply_along_axis(np.linalg.norm, axis, array)


# this is faster than using flat array
@njit(parallel=True)
def _random_normal(
    m: np.ndarray,
    v: np.ndarray,
    n_samples: int = 1,
):
    """
    Sample number from normal distribution.

    Parameters
    ----------
    m
        Mean vector.
    v
        Variance vector.
    n_samples
        Number of samples to be generated

    Returns
    -------
    :class:`numpy.ndarray`
        `(n_samples x m.shape[0])` array from normal distribution.
    """

    assert m.ndim == 1, "Means are not 1 dimensional."
    assert m.shape == v.shape, "Means and variances have different shape."

    if n_samples == 1:
        return np.expand_dims(
            np.array([np.random.normal(m[i], v[i]) for i in prange(m.shape[0])]), 0
        )

    return np.array(
        [
            [np.random.normal(m[i], v[i]) for _ in prange(n_samples)]
            for i in prange(m.shape[0])
        ]
    ).T


@njit(**jit_kwargs)
def _get_probs_for_zero_vec(size: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get a vector with uniform probability and a vector of zeros.

    Parameters
    ----------
    size
        Size of the vector.

    Returns
    -------
    :class:`numpy.ndarray`, :class:`numpy.ndarray`
        The probability and variance vectors.
    """

    # float32 doesn't have enough precision
    return (
        np.ones(size, dtype=np.float64) / size,
        np.zeros(size, dtype=np.float64),
    )


@njit(**jit_kwargs)
def _predict_transition_probabilities_numpy(
    X: np.ndarray, W: np.ndarray, softmax_scale: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute a categorical distribution based on correlation between rows in ``W`` and ``X``.

    We usually identify ``X`` with a velocity vector and ``W`` as the matrix storing transcriptomic
    displacements of the current reference cell to its nearest neighbors. For the backward process, ``X``
    is a matrix as well, storing the velocity vectors of all nearest neighbors.

    Parameters
    ----------
    X
        Either vector of shape `(n_features,)` or matrix of shape `(n_samples x n_features)`.
    W
        Weight matrix of shape `(n_samples x n_features)`.
    softmax_scale
        Scaling factor for softmax activation function.

    Returns
    --------
    :class:`scipy.sparse.csr_matrix`, :class:`scipy.sparse.csr_matrix`
        The probability and pearson correlation matrices.
    """

    # mean centering + cosine correlation
    W -= np.expand_dims(np_mean(W, axis=1), axis=1)
    W_norm = norm(W, axis=1)

    if X.shape[0] == 1:
        X = X - np.mean(X)
        X_norm = np.linalg.norm(X)

        denom = X_norm * W_norm
        mask = denom == 0
        denom[mask] = 1

        # pearson correlation
        x = W.dot(X[0]) / denom
    else:
        assert X.shape[0] == W.shape[0]
        X = X - np.expand_dims(np_mean(X, axis=1), axis=1)
        X_norm = norm(X, axis=1)

        denom = X_norm * W_norm
        mask = denom == 0
        denom[mask] = 1

        # pearson correlation
        x = np.array([np.dot(X[i], W[i]) for i in range(X.shape[0])]) / denom

    numerator = x * softmax_scale
    numerator = np.exp(numerator - np.nanmax(numerator))
    numerator = np.where(mask, 0, numerator)  # essential

    return numerator / np.nansum(numerator), x


def _filter_kwargs(fn: Callable, **kwargs) -> dict:
    """
    Filter keyword arguments.

    Parameters
    ----------
    fn
        Function for which to filter keyword arguments.
    **kwargs
        Keyword arguments to filter

    Returns
    -------
    dict
        Filtered keyword arguments for the given function.
    """

    sig = signature(fn).parameters
    return {k: v for k, v in kwargs.items() if k in sig}


def _reconstruct_one(
    data: np.ndarray,
    mat: csr_matrix,
    ixs: Optional[np.ndarray] = None,
) -> Tuple[csr_matrix, csr_matrix]:
    """
    Transform :class:`numpy.ndarray` into :class:`scipy.sparse.csr_matrix`.

    Parameters
    ----------
    data
        Array of shape `(2 x number_of_nnz)`.
    mat
        The original sparse matrix.
    ixs
        Indices that were used to sort the data.

    Returns
    -------
    :class:`scipy.sparse.csr_matrix`, :class:`scipy.sparse.csr_matrix`
        The probability and correlation matrix.
    """

    assert data.ndim == 2 and data.shape == (
        2,
        mat.nnz,
    ), f"Dimension or shape mismatch: `{data.shape}`, `{2, mat.nnz}`."

    aixs = None
    if ixs is not None:
        aixs = np.argsort(ixs)
        assert len(ixs) == mat.shape[0], f"Shape mismatch: `{ixs.shape}`, `{mat.shape}`"
        mat = mat[ixs]

    # strange bug happens when no copying and eliminating zeros from cors (it's no longer row-stochastic)
    # only happens when using numba
    probs = csr_matrix((data[0], mat.indices, mat.indptr))
    cors = csr_matrix((data[1], mat.indices, mat.indptr))

    if aixs is not None:
        assert (
            len(aixs) == probs.shape[0]
        ), f"Shape mismatch: `{ixs.shape}`, `{probs.shape}`."
        probs, cors = probs[aixs], cors[aixs]

    probs.eliminate_zeros()
    cors.eliminate_zeros()

    close_to_1 = np.isclose(probs.sum(1), 1.0)
    if not np.all(close_to_1):
        raise ValueError(
            f"Matrix is not row-stochastic. "
            f"The following rows don't sum to 1: `{list(np.where(~close_to_1)[0])}`."
        )

    return probs, cors


@njit(**jit_kwargs)
def _calculate_starts(indptr: np.ndarray, ixs) -> np.ndarray:
    """
    Get the position where to put the data.

    Parameters
    ----------
    indptr
        Pointer of indices from :class:`scipy.sparse.csr_matrix`.
    ixs
        Row indices for which to calculate the starts.

    Returns
    -------
    :class:`numpy.ndarray`
        The starting positions.
    """

    starts = np.cumsum(indptr[ixs + 1] - indptr[ixs])
    return np.hstack((np.array([0], dtype=starts.dtype), starts))
