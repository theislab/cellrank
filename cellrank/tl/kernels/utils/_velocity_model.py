from typing import Any, Tuple, Union, Callable, Optional

from abc import ABC, abstractmethod
from enum import auto
from math import fsum

from cellrank.tl._enum import _DEFAULT_BACKEND, ModeEnum, Backend_t
from cellrank.ul._parallelize import parallelize
from cellrank.tl.kernels._utils import _random_normal, _calculate_starts
from cellrank.tl.kernels._base_kernel import _RTOL
from cellrank.tl.kernels.utils._similarity_scheme import Similarity

import numpy as np
from numba import prange
from scipy.sparse import spmatrix, csr_matrix

__all__ = ("Deterministic", "Stochastic", "MonteCarlo")


class VelocityModel(ModeEnum):
    DETERMINISTIC = auto()
    STOCHASTIC = auto()
    MONTE_CARLO = auto()


class BackwardMode(ModeEnum):
    TRANSPOSE = auto()
    NEGATE = auto()


class ModelABC(ABC):
    def __init__(
        self,
        conn: spmatrix,
        x: np.ndarray,
        v: np.ndarray,
        similarity: Union[Similarity, Callable[[np.ndarray], np.ndarray]],
        backward_mode: Optional[BackwardMode] = None,
        softmax_scale: float = 1.0,
        dtype: np.dtype = np.float64,
    ):
        self._conn = conn.astype(dtype, copy=False)
        self._x = x.astype(dtype, copy=False)
        self._v = v.astype(dtype, copy=False)
        self._similarity = similarity
        self._backward_mode = backward_mode
        self._softmax_scale = softmax_scale
        self._dtype = dtype

    @abstractmethod
    def _compute(self, ix: int, neighs_ixs: np.ndarray, **kwargs: Any) -> np.ndarray:
        pass

    def __call__(
        self,
        n_jobs: Optional[int] = None,
        backend: Backend_t = _DEFAULT_BACKEND,
        **kwargs: Any,
    ) -> Tuple[np.ndarray, np.ndarray]:
        ixs = self._ixs
        return parallelize(
            self._compute_helper,
            ixs,
            n_jobs=n_jobs,
            backend=backend,
            as_array=False,
            extractor=lambda data: self._reconstruct_output(
                np.concatenate(data, axis=-1), ixs
            ),
            unit=self._unit,
        )(**kwargs)

    def _compute_helper(
        self, ixs: np.ndarray, queue=None, **kwargs
    ) -> Tuple[np.ndarray, np.ndarray]:
        indptr, indices = self._conn.indptr, self._conn.indices
        starts = _calculate_starts(indptr, ixs)
        probs_logits = np.empty((2, starts[-1]), dtype=np.float64)

        for i, ix in enumerate(ixs):
            start, end = indptr[ix], indptr[ix + 1]
            neigh_ixs = indices[start:end]
            n_neigh = len(neigh_ixs)

            tmp = self._compute(ix, neigh_ixs, **kwargs)
            if np.shape(tmp) != (2, n_neigh):
                raise ValueError(
                    f"Expected row of shape `{(2, n_neigh)}`, found `{np.shape(tmp)}`."
                )

            probs_logits[:, starts[i] : starts[i] + n_neigh] = tmp
            if queue is not None:
                queue.put(1)

        if queue is not None:
            queue.put(None)

        return probs_logits

    def _reconstruct_output(
        self,
        data: np.ndarray,
        ixs: Optional[np.ndarray] = None,
    ) -> csr_matrix:
        """
        Transform :class:`numpy.ndarray` into :class:`scipy.sparse.csr_matrix`.

        Parameters
        ----------
        data
            Array of shape ``(2, nnz)``.
        ixs
            Indices that were used to sort the data.

        Returns
        -------
        TODO.
        """

        def reconstruct(data: np.ndarray) -> Tuple[csr_matrix, csr_matrix]:
            data = csr_matrix(
                (
                    np.array(data, copy=True),
                    np.array(conn.indices, copy=True),
                    np.array(conn.indptr, copy=True),
                )
            )
            if aixs is not None:
                data = data[aixs]
            data.eliminate_zeros()
            return data

        if data.shape != (2, self._conn.nnz):
            raise ValueError(
                f"Dimension or shape mismatch: `{data.shape}`, `{(2, self._conn.nnz)}`."
            )

        if ixs is None:
            aixs, conn = None, self._conn
        else:
            aixs = np.argsort(ixs)
            conn = self._conn[ixs]

        # strange bug happens when no copying and eliminating zeros from cors (it's no longer row-stochastic)
        # only happens when using numba
        probs = reconstruct(data[0])
        logits = reconstruct(data[1])
        if not np.allclose(probs.sum(1), 1.0, rtol=_RTOL):
            # TODO: row ixs?
            raise ValueError(f"Matrix is not row-stochastic.")

        return probs, logits

    @property
    def _ixs(self) -> np.ndarray:
        """Order in which to process rows."""
        ixs = np.arange(self._conn.shape[0])
        np.random.shuffle(ixs)
        return ixs

    def _uniform(self, n: int) -> Tuple[np.ndarray, np.ndarray]:
        """Uninformative probability vector when velocity is 0."""
        return np.ones(n, dtype=self._dtype) / n, np.zeros((n,), dtype=self._dtype)

    @property
    def _unit(self) -> str:
        """Progress bar unit."""
        return "cell"


class Deterministic(ModelABC):
    def _compute(self, ix: int, neigh_ixs: np.ndarray, **_: Any) -> np.ndarray:
        W = self._x[neigh_ixs, :] - self._v[ix, :]

        if self._backward_mode not in (None, BackwardMode.NEGATE):
            return self._similarity(self._v[neigh_ixs, :], -1 * W, self._softmax_scale)

        v = self._v[ix]
        if np.all(v == 0):
            return self._uniform(len(neigh_ixs))

        if self._backward_mode == BackwardMode.NEGATE:
            v *= -1.0

        return self._similarity(v[None, :], W, self._softmax_scale)


class Stochastic(ModelABC):
    def __init__(
        self,
        conn: spmatrix,
        x: np.ndarray,
        exp: np.ndarray,
        var: np.ndarray,
        dtype: np.dtype = np.float64,
        **kwargs: Any,
    ):
        super().__init__(conn, x, exp, **kwargs)
        if not hasattr(self._similarity, "hessian"):
            raise AttributeError("Similarity scheme doesn't have a `hessian` function.")
        self._var = var.astype(dtype, copy=False)

    def _compute(
        self, ix: int, nbhs_ixs: np.ndarray, **_: Any
    ) -> Tuple[np.ndarray, np.ndarray]:
        v = self._v[ix]
        n_neigh, n_feat = len(nbhs_ixs), self._x.shape[1]

        if np.all(v == 0):
            return self._uniform(n_neigh)
        # compute the Hessian tensor, and turn it into a matrix that has the diagonal elements in its rows
        W = self._x[nbhs_ixs, :] - self._x[ix, :]

        H = self._similarity.hessian(v, W, self._softmax_scale)
        if H.shape == (n_neigh, n_feat):
            H_diag = H
        elif H.shape == (n_neigh, n_feat, n_feat):
            H_diag = np.array([np.diag(h) for h in H])
        else:
            raise ValueError(
                f"Expected full Hessian matrix of shape `{(n_neigh, n_feat, n_feat)}` "
                f"or its diagonal of shape `{(n_neigh, n_feat)}`, found `{H.shape}`."
            )

        # compute zero order term
        p_0, c = self._similarity(v[None, :], W, self._softmax_scale)

        # compute second order term (note that the first order term cancels)
        p_2 = 0.5 * H_diag.dot(self._var[ix])

        # combine both to give the second order Taylor approximation.
        # Can sometimes be negative because we neglected higher order terms, so force it to be non-negative
        p = np.clip(p_0 + p_2, a_min=0, a_max=1)

        nan_mask = np.isnan(p)
        p[nan_mask] = 0
        c[nan_mask] = 0

        if np.all(p == 0):
            return self._uniform(n_neigh)

        sum_ = fsum(p)
        if not np.isclose(sum_, 1.0, rtol=_RTOL):
            p[~nan_mask] = p[~nan_mask] / sum_

        return p, c

    @property
    def _ixs(self) -> np.ndarray:
        # less recompilation for JAX
        return np.argsort((self._conn != 0).sum(1).A1)[::-1]


class MonteCarlo(ModelABC):
    def __init__(
        self,
        conn: spmatrix,
        x: np.ndarray,
        exp: np.ndarray,
        var: np.ndarray,
        n_samples: int = 1,
        seed: Optional[int] = None,
        dtype: np.dtype = np.float64,
        **kwargs: Any,
    ):
        super().__init__(conn, x, exp, **kwargs)
        self._var = var.astype(dtype, np.float64)
        self._n_samples = n_samples
        np.random.seed(seed)

    def _compute(
        self, ix: int, nbhs_ixs: np.ndarray, **_: Any
    ) -> Tuple[np.ndarray, np.ndarray]:
        # fmt: off
        n_neigh = len(nbhs_ixs)
        W = self._x[nbhs_ixs, :] - self._x[ix, :]

        samples = _random_normal(self._x[ix], self._var[ix], n_samples=self._n_samples)

        probs = np.zeros((n_neigh,), dtype=np.float64)
        logits = np.zeros((n_neigh,), dtype=np.float64)
        for j in prange(self._n_samples):
            p, l = self._similarity(np.atleast_2d(samples[j]), W, self._softmax_scale)
            probs += p
            logits += l

        return probs / self._n_samples, logits / self._n_samples
        # fmt: on

    def _unit(self) -> str:
        return "sample"
