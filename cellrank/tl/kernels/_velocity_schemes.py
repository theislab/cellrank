from abc import ABC, abstractmethod
from typing import Any, Tuple
from functools import partial

import numpy as np
from numba import njit

from cellrank.ul._docs import d
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


class Hessian(ABC):  # noqa: D101
    @d.get_full_description(base="hessian")
    @d.get_sections(base="hessian", sections=["Parameters", "Returns"])
    @abstractmethod
    def hessian(
        self, v: np.ndarray, D: np.ndarray, softmax_scale: float = 1.0
    ) -> np.ndarray:
        """
        Compute the Hessian.

        Parameters
        ----------
        v
            Array of shape ``(n_genes,)`` containing the velocity vector.
        D
            Array of shape ``(n_neighbors, n_genes)`` corresponding to the transcriptomic displacement of the current
            cell with respect to ist nearest neighbors.
        softmax_scale
            Scaling factor for the softmax function.

        Returns
        -------
        :class:`numpy.ndarray`
            The full Hessian of shape ``(n_neighbors, n_genes, n_genes)`` or only its diagonal of shape
            ``(n_neighbors, n_genes)``.

        Developer notes
        ---------------
        This class should be used in conjunction with any class that implements the ``__call__`` method and should
        compute the Hessian of such function. This it does not need to handle the case of a backward process,
        i.e. when the velocity vector `v` is of shape ``(n_genes, n_neighbors)``.

        If using :mod:`jax` to compute the Hessian, please specify a class attribute ``__use_jax__ = True``.
        """


class SimilaritySchemeABC(ABC):
    """Base class for all similarity schemes."""

    @d.get_full_description(base="sim_scheme")
    @d.get_sections(base="sim_scheme", sections=["Parameters", "Returns"])
    @abstractmethod
    def __call__(
        self, v: np.ndarray, D: np.ndarray, softmax_scale: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute transition probability of a cell to its nearest neighbors using RNA velocity.

        Parameters
        ----------
        v
            Array of shape ``(n_genes,)`` or ``(n_neighbors, n_genes)`` containing the velocity vector(s).
            The second case is used for the backward process.
        D
            Array of shape ``(n_neighbors, n_genes)`` corresponding to the transcriptomic displacement of the current
            cell with respect to ist nearest neighbors.
        softmax_scale
            Scaling factor for the softmax function.

        Returns
        -------
        :class:`numpy.ndarray`, :class:`numpy.ndarray`
            The probability and logits arrays of shape ``(n_neighbors,)``.
        """

    def __repr__(self):
        return f"<{self.__class__.__name__}>"

    def __str__(self):
        return repr(self)


class SimilaritySchemeHessian(SimilaritySchemeABC, Hessian):
    """
    Base class for all similarity schemes as defined in [Li2020]_.

    Parameters
    ----------
    center_mean
        Whether to center the velocity vector(s) and the transcriptomic displacement matrix.
    scale_by_norm
        Whether to scale the velocity vector(s) and the transcriptomic displacement matrix by their norms.
    """

    __use_jax__: bool = True

    def __init__(self, center_mean: bool, scale_by_norm: bool):
        self._center_mean = center_mean
        self._scale_by_norm = scale_by_norm

    @d.dedent
    def __call__(
        self, v: np.ndarray, D: np.ndarray, softmax_scale: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        %(sim_scheme.full_desc)s

        Parameters
        ----------
        %(sim_scheme.parameters)s

        Returns
        -------
        %(sim_scheme.returns)s
        """  # noqa: D400
        return _predict_transition_probabilities_numpy(
            v, D, softmax_scale, self._center_mean, self._scale_by_norm
        )

    @d.dedent
    def hessian(
        self, v: np.ndarray, D: np.ndarray, softmax_scale: float = 1.0
    ) -> np.ndarray:
        """
        %(hessian.full_desc)s

        Parameters
        ----------
        %(hessian.parameters)s

        Returns
        -------
        %(hessian.returns)s
        """  # noqa: D400
        return _predict_transition_probabilities_jax_H(
            v, D, softmax_scale, self._center_mean, self._scale_by_norm
        )

    def __getattribute__(self, name: str) -> Any:
        if name == "hessian" and self.__use_jax__ and not _HAS_JAX:
            raise NotImplementedError("No `jax` installation found.")
        return super().__getattribute__(name)


class DotProductScheme(SimilaritySchemeHessian):
    r"""
    Dot product scheme as defined in eq. (4.9) of [Li2020]_.

        :math:`v(s_i, s_j) = g(\delta_{i, j}^T v_i)`

    where :math:`v_i` is the velocity vector of cell :math:`i`, :math:`\delta_{i, j}` corresponds to the transcriptional
    displacement between cells :math:`i` and :math:`j` and :math:`g` is a softmax function with some scaling parameter.
    """

    def __init__(self):
        super().__init__(center_mean=False, scale_by_norm=False)


class CosineScheme(SimilaritySchemeHessian):
    r"""
    Cosine similarity scheme as defined in eq. (4.7) of [Li2020]_.

        :math:`v(s_i, s_j) = g(cos(\delta_{i, j}, v_i))`

    where :math:`v_i` is the velocity vector of cell :math:`i`, :math:`\delta_{i, j}` corresponds to the transcriptional
    displacement between cells :math:`i` and :math:`j` and :math:`g` is a softmax function with some scaling parameter.
    """

    def __init__(self):
        super().__init__(center_mean=False, scale_by_norm=True)


class CorrelationScheme(SimilaritySchemeHessian):
    r"""
    Pearson correlation scheme as defined in eq. (4.8) of [Li2020]_.

        :math:`v(s_i, s_j) = g(corr(\delta_{i, j}, v_i))`

    where :math:`v_i` is the velocity vector of cell :math:`i`, :math:`\delta_{i, j}` corresponds to the transcriptional
    displacement between cells :math:`i` and :math:`j` and :math:`g` is a softmax function with some scaling parameter.
    """

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
