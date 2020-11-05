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

    @partial(jit, static_argnums=(3,))
    def _predict_transition_probabilities_jax(
        X: np.ndarray,
        W: np.ndarray,
        softmax_scale: float = 1.0,
        center_mean: bool = True,
    ):
        if center_mean:
            # pearson correlation, otherwise cosine
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

        return numerator / jnp.nansum(numerator)

    _predict_transition_probabilities_jax_H = hessian(
        _predict_transition_probabilities_jax
    )

    _HAS_JAX = True
except ImportError:
    _HAS_JAX = False

    hessian = jit = lambda _: _

    def _predict_transition_probabilities_jax(*args, **kwargs):
        raise NotImplementedError("No `jax installation found`.")

    _predict_transition_probabilities_jax_H = _predict_transition_probabilities_jax


@njit(**jit_kwargs)
def _predict_transition_probabilities_numpy(
    X: np.ndarray,
    W: np.ndarray,
    softmax_scale: float = 1.0,
    center_mean: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    if center_mean:
        # pearson correlation, otherwise cosine
        W -= np.expand_dims(np_mean(W, axis=1), axis=1)
    W_norm = norm(W, axis=1)

    if X.shape[0] == 1:
        if center_mean:
            X = X - np.mean(X)
        X_norm = np.linalg.norm(X)

        denom = X_norm * W_norm
        mask = denom == 0
        denom[mask] = 1

        # pearson correlation
        x = W.dot(X[0]) / denom
    else:
        assert X.shape[0] == W.shape[0], "Wrong shape."
        if center_mean:
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


class SimilaritySchemeMeta(ABCMeta):  # noqa: D101
    def __new__(cls, clsname, superclasses, attributedict):  # noqa: D102
        if (
            attributedict.pop("__use_jax__", False)
            and not _HAS_JAX
            and "hessian" in attributedict
        ):
            del attributedict["hessian"]
        res = super().__new__(cls, clsname, superclasses, attributedict)

        return res


class SimilaritySchemeABC(ABC, metaclass=SimilaritySchemeMeta):  # noqa: D101
    __use_jax__: bool = True

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
        pass

    @abstractmethod
    def _hessian(
        self, X: np.ndarray, W: np.ndarray, softmax_scale: float = 1.0
    ) -> np.ndarray:
        pass

    def hessian(
        self, X: np.ndarray, W: np.ndarray, softmax_scale: float = 1.0
    ) -> np.ndarray:
        """TODO."""
        return self._hessian(X, W, softmax_scale)

    def __repr__(self):
        return f"<{self.__class__.__name__}>"

    def __str__(self):
        return repr(self)


# TODO: rename me
class Scheme_(SimilaritySchemeABC):  # noqa: D101
    # TODO: refrain from saving state inside
    def __init__(self, center_mean: bool):
        self._center_mean = center_mean

    def __call__(
        self, X: np.ndarray, W: np.ndarray, softmax_scale: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """TODO: docrep."""
        return _predict_transition_probabilities_numpy(
            X, W, softmax_scale, self._center_mean
        )

    def _hessian(
        self, X: np.ndarray, W: np.ndarray, softmax_scale: float = 1.0
    ) -> np.ndarray:
        return _predict_transition_probabilities_jax_H(
            X, W, softmax_scale, self._center_mean
        )


class CosineScheme(Scheme_):  # noqa: D101 TODO
    def __init__(self):
        super().__init__(center_mean=False)


class CorrelationScheme(Scheme_):  # noqa: D101 TODO
    def __init__(self):
        super().__init__(center_mean=True)


class DotProductScheme(SimilaritySchemeABC):  # noqa: D101 TODO
    # TODO: refrain from saving state inside
    def __init__(self):
        def _hessian_impl(
            X: np.ndarray, W: np.ndarray, softmax_scale=1.0
        ) -> np.ndarray:
            numerator = W.dot(X) * softmax_scale
            numerator = jnp.exp(numerator - jnp.nanmax(numerator))

            return numerator / jnp.nansum(numerator)

        # TODO: add some activation function g?
        self._hessian_impl = hessian(jit(_hessian_impl))

    # hardly any point in jitting this
    def __call__(
        self, X: np.ndarray, W: np.ndarray, softmax_scale: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """TODO: docrep."""
        if X.shape[0] == 1:
            x = W.dot(X[0])
        else:
            x = np.array([np.dot(X[i], W[i]) for i in range(X.shape[0])])

        numerator = x * softmax_scale
        numerator = np.exp(numerator - np.nanmax(numerator))

        return numerator / np.nansum(numerator), x

    def _hessian(
        self, X: np.ndarray, W: np.ndarray, softmax_scale: float = 1.0
    ) -> np.ndarray:
        return self._hessian_impl(X, W, softmax_scale)


@valuedispatch
def _create_scheme(scheme: Scheme, *_args, **_kwargs) -> SimilaritySchemeABC:
    raise NotImplementedError(scheme)


@_create_scheme.register(Scheme.DOT_PRODUCT)
def _():
    return DotProductScheme()


@_create_scheme.register(Scheme.COSINE)
def _():
    return CosineScheme()


@_create_scheme.register(Scheme.CORRELATION)
def _():
    return CorrelationScheme()
