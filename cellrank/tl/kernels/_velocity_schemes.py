# -*- coding: utf-8 -*-
from abc import ABC, ABCMeta, abstractmethod
from typing import Tuple
from functools import partial

import numpy as np
from numba import njit

from cellrank.ul._utils import valuedispatch
from cellrank.tl._constants import ModeEnum
from cellrank.tl.kernels._utils import norm, np_mean, jit_kwargs


class Scheme(ModeEnum):  # noqa: D101
    DOT_PRODUCT = "dot_product"
    COSINE = "cosine"
    CORRELATION = "correlation"


try:
    import jax.numpy as jnp
    from jax import jit, hessian

    @jit
    def _softmax_jax(x: np.ndarray, softmax_scale) -> np.ndarray:
        numerator = x * softmax_scale
        numerator = jnp.exp(numerator - jnp.max(numerator))
        return numerator / jnp.sum(numerator)

    @jit
    def _softmax_masked_jax(
        x: np.ndarray, mask: np.ndarray, softmax_scale
    ) -> np.ndarray:
        numerator = x * softmax_scale
        numerator = jnp.exp(numerator - jnp.nanmax(numerator))
        numerator = jnp.where(mask, 0, numerator)  # essential

        return numerator / jnp.nansum(numerator)

    @partial(jit, static_argnums=(3, 4))
    def _predict_transition_probabilities_jax(
        X: np.ndarray,
        W: np.ndarray,
        softmax_scale: float = 1.0,
        center_mean: bool = True,
        scale_by_norm: bool = True,
    ):
        if center_mean:
            # pearson correlation, otherwise cosine
            W -= W.mean(axis=1)[:, None]
            X -= X.mean()

        if scale_by_norm:
            denom = jnp.linalg.norm(X) * jnp.linalg.norm(W, axis=1)
            mask = jnp.isclose(denom, 0)
            denom = jnp.where(jnp.isclose(denom, 0), 1, denom)  # essential
            return _softmax_masked_jax(W.dot(X) / denom, mask, softmax_scale)

        return _softmax_jax(W.dot(X), softmax_scale)

    _predict_transition_probabilities_jax_H = hessian(
        _predict_transition_probabilities_jax
    )

    _HAS_JAX = True
except ImportError:
    _HAS_JAX = False

    hessian = jit = lambda _: _

    def _predict_transition_probabilities_jax(*_args, **_kwargs):
        raise NotImplementedError("No `jax` installation found.")

    _predict_transition_probabilities_jax_H = _predict_transition_probabilities_jax


@njit(**jit_kwargs)
def _softmax(x: np.ndarray, softmax_scale: float) -> Tuple[np.ndarray, np.ndarray]:
    numerator = x * softmax_scale
    numerator = np.exp(numerator - np.max(numerator))
    return numerator / np.sum(numerator), x


@njit(**jit_kwargs)
def _softmax_masked(
    x: np.ndarray, mask: np.ndarray, softmax_scale: float
) -> Tuple[np.ndarray, np.ndarray]:
    numerator = x * softmax_scale
    numerator = np.exp(numerator - np.nanmax(numerator))
    numerator = np.where(mask, 0, numerator)  # essential

    return numerator / np.nansum(numerator), x


@njit(parallel=False, **jit_kwargs)
def _predict_transition_probabilities_numpy(
    X: np.ndarray,
    W: np.ndarray,
    softmax_scale: float = 1.0,
    center_mean: bool = True,
    scale_by_norm: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    if center_mean:
        # pearson correlation
        W -= np.expand_dims(np_mean(W, axis=1), axis=1)

    if X.shape[0] == 1:
        if center_mean:
            # pearson correlation
            X = X - np.mean(X)

        if scale_by_norm:
            # cosine or pearson correlation
            denom = np.linalg.norm(X) * norm(W, axis=1)
            mask = denom == 0
            denom[mask] = 1
            return _softmax_masked(W.dot(X[0]) / denom, mask, softmax_scale)

        return _softmax(W.dot(X[0]), softmax_scale)

    assert X.shape[0] == W.shape[0], "Wrong shape."

    if center_mean:
        X = X - np.expand_dims(np_mean(X, axis=1), axis=1)

    if scale_by_norm:
        denom = norm(X, axis=1) * norm(W, axis=1)
        mask = denom == 0
        denom[mask] = 1
        return _softmax_masked(
            np.array([np.dot(X[i], W[i]) for i in range(X.shape[0])]) / denom,
            mask,
            softmax_scale,
        )

    return _softmax(
        np.array([np.dot(X[i], W[i]) for i in range(X.shape[0])]), softmax_scale
    )


class SimilaritySchemeMeta(ABCMeta):  # noqa: D101
    def __new__(cls, clsname, superclasses, attributedict):  # noqa: D102
        res = super().__new__(cls, clsname, superclasses, attributedict)
        if (
            attributedict.pop("__use_jax__", False)
            and not _HAS_JAX
            and "hessian" in attributedict
        ):
            delattr(res, "hessian")

        return res


class Hessian(ABC, metaclass=SimilaritySchemeMeta):  # noqa: D101
    @abstractmethod
    def hessian(
        self, X: np.ndarray, W: np.ndarray, softmax_scale: float = 1.0
    ) -> np.ndarray:
        """TODO."""
        pass


class SimilaritySchemeABC(ABC):
    """TODO."""

    @abstractmethod
    def __call__(
        self, X: np.ndarray, W: np.ndarray, softmax_scale: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        # TODO: update docs
        """
        Compute a categorical distribution based on correlation between rows in ``W`` and ``X``.

        We usually identify ``X`` with a velocity vector and ``W`` as the matrix storing transcriptomic
        displacements of the current reference cell to its nearest neighbors. For the backward process, ``X``
        is a matrix as well, storing the velocity vectors of all nearest neighbors.

        Parameters
        ----------
        X
            Either vector of shape `(n_features,)` or matrix of shape ``(n_samples, n_features)``.
        W
            Weight matrix of shape ``(n_samples, n_features)``.
        softmax_scale
            Scaling factor for softmax activation function.

        Returns
        --------
        :class:`scipy.sparse.csr_matrix`, :class:`scipy.sparse.csr_matrix`
            The probability and pearson correlation matrices.
        """
        pass

    def __repr__(self):
        return f"<{self.__class__.__name__}>"

    def __str__(self):
        return repr(self)


class SimilarityScheme(SimilaritySchemeABC, Hessian):  # noqa: D101
    __use_jax__: bool = True

    def __init__(self, center_mean: bool, scale_by_norm: bool):
        self._center_mean = center_mean
        self._scale_by_norm = scale_by_norm

    def __call__(
        self, X: np.ndarray, W: np.ndarray, softmax_scale: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """TODO."""
        return _predict_transition_probabilities_numpy(
            X, W, softmax_scale, self._center_mean, self._scale_by_norm
        )

    def hessian(
        self, X: np.ndarray, W: np.ndarray, softmax_scale: float = 1.0
    ) -> np.ndarray:
        """TODO."""
        return _predict_transition_probabilities_jax_H(
            X, W, softmax_scale, self._center_mean, self._scale_by_norm
        )


class DotProductScheme(SimilarityScheme):  # noqa: D101
    def __init__(self):
        super().__init__(center_mean=False, scale_by_norm=False)


class CosineScheme(SimilarityScheme):  # noqa: D101 TODO
    def __init__(self):
        super().__init__(center_mean=False, scale_by_norm=True)


class CorrelationScheme(SimilarityScheme):  # noqa: D101 TODO
    def __init__(self):
        super().__init__(center_mean=True, scale_by_norm=True)


@valuedispatch
def _get_scheme(scheme: Scheme, *_args, **_kwargs) -> SimilaritySchemeABC:
    raise NotImplementedError(scheme)


@_get_scheme.register(Scheme.DOT_PRODUCT)
def _():
    return DotProductScheme()


@_get_scheme.register(Scheme.COSINE)
def _():
    return CosineScheme()


@_get_scheme.register(Scheme.CORRELATION)
def _():
    return CorrelationScheme()
