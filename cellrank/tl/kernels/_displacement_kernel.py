from typing import Any, Tuple, Union, Callable, Optional
from typing_extensions import Literal

from abc import ABC, abstractmethod
from copy import copy
from enum import auto
from math import fsum

from anndata import AnnData
from cellrank import logging as logg
from cellrank.tl._enum import _DEFAULT_BACKEND, ModeEnum
from cellrank.ul._docs import d, inject_docs
from cellrank.ul._utils import valuedispatch
from cellrank.tl.kernels._bk import BidirectionalKernel
from cellrank.ul._parallelize import parallelize
from cellrank.tl.kernels._utils import (
    prange,
    np_mean,
    _filter_kwargs,
    _random_normal,
    _reconstruct_one,
    _calculate_starts,
    _get_probs_for_zero_vec,
)
from cellrank.tl.kernels._mixins import ConnectivityMixin
from scvelo.preprocessing.moments import get_moments
from cellrank.tl.kernels._base_kernel import _RTOL
from cellrank.tl.kernels._velocity_schemes import Scheme, SimilaritySchemeABC

import numpy as np
from scipy.sparse import issparse, spmatrix, csr_matrix


class VelocityMode(ModeEnum):  # noqa: D101
    DETERMINISTIC = auto()
    STOCHASTIC = auto()
    SAMPLING = auto()
    MONTE_CARLO = auto()


class BackwardMode(ModeEnum):  # noqa: D101
    TRANSPOSE = auto()
    NEGATE = auto()


@d.dedent
class DisplacementKernel(ConnectivityMixin, BidirectionalKernel):
    def __init__(
        self,
        adata: AnnData,
        key1: Optional[str],
        key2: Optional[str],
        key1_displacement: bool = True,
        key2_displacement: bool = True,
        backward: bool = False,
        check_connectivity: bool = False,
        **kwargs: Any,
    ):
        super().__init__(
            adata,
            key1=key1,
            key2=key2,
            key1_displacement=key1_displacement,
            key2_displacement=key2_displacement,
            backward=backward,
            check_connectivity=check_connectivity,
            **kwargs,
        )
        self._vkey = key2  # for copy
        self._xkey = key1

    def _extract_layer(
        self, key: Optional[str] = None, dtype: np.dtype = np.float64
    ) -> np.ndarray:
        if key in (None, "X"):
            data = self.adata.X
        elif key in self.adata.layers:
            data = self.adata.layers[key]
        else:
            raise KeyError("TODO.")

        data = data.astype(dtype)
        return data.toarray() if issparse(data) else data

    def _read_from_adata(
        self,
        key1: str,
        key2: str,
        key1_displacement: bool = True,
        key2_displacement: bool = True,
        **kwargs: Any,
    ) -> None:
        super()._read_from_adata(**kwargs)

        self._data1 = self._extract_layer(key1)
        self._data2 = self._extract_layer(key2)
        self._displacement1 = key1_displacement
        self._displacement2 = key2_displacement

    def _get_scheme_and_model(
        self,
        mode: Union[str, VelocityMode],
        scheme,
        n_samples: int = 1000,
        backend: str = _DEFAULT_BACKEND,
    ) -> Tuple[VelocityMode, Callable, str]:
        mode = VelocityMode(mode)
        if isinstance(scheme, str):
            scheme = SimilaritySchemeABC.create(Scheme(scheme))
        elif not callable(scheme):
            raise TypeError(
                f"Expected `scheme` to be a function, found `{type(scheme).__name__}`."
            )

        if self.backward and mode != VelocityMode.DETERMINISTIC:
            logg.warning(
                f"Mode `{mode!r}` is currently not supported for backward process. "
                f"Using to `mode={VelocityMode.DETERMINISTIC!r}`"
            )
            mode = VelocityMode.DETERMINISTIC

        if mode == VelocityMode.STOCHASTIC and not hasattr(scheme, "hessian"):
            mode = VelocityMode.MONTE_CARLO
            logg.warning(
                f"Unable to detect a method for Hessian computation. If using one of the "
                f"predefined similarity functions, consider installing `jax` as "
                f"`pip install jax jaxlib`. Using `mode={mode!r}` and `n_samples={n_samples}`"
            )

        if mode == VelocityMode.MONTE_CARLO and n_samples == 1:
            mode = VelocityMode.SAMPLING

        if mode != VelocityMode.STOCHASTIC and backend == "multiprocessing":
            # TODO(michalk8): this is because on jitting and pickling (cloudpickle, used by loky, handles it correctly)
            backend = _DEFAULT_BACKEND
            logg.warning(
                f"Multiprocessing backend is supported only for `mode={VelocityMode.STOCHASTIC!r}`. Using `{backend!r}`"
            )

        return mode, scheme, backend

    def _estimate_softmax_scale(self, **kwargs) -> float:
        1.0 / np.median(np.abs(cmat.data))
        return 0

    @inject_docs(m=VelocityMode, b=BackwardMode, s=Scheme)  # don't swap the order
    @d.dedent
    def compute_transition_matrix(
        self,
        mode: Literal[
            "deterministic", "stochastic", "sampling", "monte_carlo"
        ] = VelocityMode.DETERMINISTIC,
        backward_mode: Literal["transpose", "negate"] = BackwardMode.TRANSPOSE,
        scheme: Union[
            Literal["dot_product", "cosine", "correlation"], Callable
        ] = Scheme.CORRELATION,
        softmax_scale: Optional[float] = None,
        n_samples: int = 1000,
        seed: Optional[int] = None,
        **kwargs: Any,
    ) -> "DisplacementKernel":
        """
        Compute transition matrix based on velocity directions on the local manifold.

        For each cell, infer transition probabilities based on the cell's velocity-extrapolated cell state and the
        cell states of its *K* nearest neighbors.

        Parameters
        ----------
        %(velocity_mode)s
        %(velocity_backward_mode)s
        %(softmax_scale)s
        %(velocity_scheme)s
        n_samples
            Number of bootstrap samples when ``mode = {m.MONTE_CARLO!r}``.
        seed
            Random seed when ``mode = {m.MONTE_CARLO!r}``.
        %(parallel)s

        Returns
        -------
        Self and updates :attr:`transition_matrix` and :attr:`params`.
        """
        mode, scheme, backend = self._get_scheme_and_model(
            mode, scheme, n_samples, backend=kwargs.pop("backend", _DEFAULT_BACKEND)
        )
        backward_mode = BackwardMode(backward_mode)

        start = logg.info(f"Computing transition matrix using `{mode!r}` mode")

        # fmt: off
        params = {"softmax_scale": softmax_scale, "mode": mode, "seed": seed, "scheme": str(scheme)}
        if self.backward:
            params["bwd_mode"] = str(backward_mode)
        if self._reuse_cache(params, time=start):
            return self
        # fmt: on

        if softmax_scale is None:
            softmax_scale = self._estimate_softmax_scale()
            params["softmax_scale"] = softmax_scale
            logg.info(f"Setting `softmax_scale={softmax_scale:.4f}`")

        if softmax_scale is None:
            logg.info(
                f"Estimating `softmax_scale` using `{VelocityMode.DETERMINISTIC!r}` mode"
            )
            _, cmat = _dispatch_computation(
                VelocityMode.DETERMINISTIC,
                scheme=scheme,
                conn=self._conn,
                expression=self._gene_expression,
                velocity=self._velocity,
                expectation=velocity_expectation,
                variance=velocity_variance,
                softmax_scale=1.0,
                backward=self.backward,
                backward_mode=backward_mode,
                n_samples=n_samples,
                seed=seed,
                backend=backend,
                **kwargs,
            )

        tmat, cmat = _dispatch_computation(
            mode,
            scheme=scheme,
            conn=self._conn,
            expression=self._gene_expression,
            velocity=self._velocity,
            expectation=velocity_expectation,
            variance=velocity_variance,
            softmax_scale=softmax_scale,
            backward=self.backward,
            backward_mode=backward_mode,
            n_samples=n_samples,
            seed=seed,
            backend=backend,
            **kwargs,
        )
        self.transition_matrix = tmat

        logg.info("    Finish", time=start)

        return self

    def __invert__(self) -> "BidirectionalMixin":
        pass

    @d.dedent
    def copy(self) -> "DisplacementKernel":
        """%(copy)s"""  # noqa
        vk = DisplacementKernel(
            self.adata,
            backward=self.backward,
            key2=self._vkey,
            key1=self._xkey,
            gene_subset=self._gene_subset,
        )
        vk._params = copy(self.params)
        vk._cond_num = self.condition_number
        vk._transition_matrix = copy(self._transition_matrix)

        return vk


class TmatCalculator(ABC):
    def __init__(
        self,
        conn: spmatrix,
        x: np.ndarray,
        v: np.ndarray,
        backward_mode: Optional[BackwardMode] = None,
        softmax_scale: float = 1.0,
    ):
        self._conn = conn
        self._x = x
        self._v = v
        self._similarity = ...
        self._backward_mode = backward_mode
        self._softmax_scale = softmax_scale

    # TODO: reconstruct one

    @abstractmethod
    def _compute_transtion_matrix(
        self, ix: int, neighs_ixs: np.ndarray
    ) -> Union[np.ndarray, spmatrix]:
        pass

    def __call__(
        self,
        *args: Any,
        n_jobs: Optional[int] = None,
        backend: str = _DEFAULT_BACKEND,
        **kwargs: Any,
    ) -> np.ndarray:
        tmat = self._compute_transtion_matrix(*args, **kwargs)
        # TODO
        return tmat


class Deterministic(TmatCalculator):
    def _compute_transtion_matrix(
        self, ix: int, neigh_ixs: np.ndarray
    ) -> Union[np.ndarray, spmatrix]:
        W = self._x[neigh_ixs, :] - self._v[ix, :]

        if self._backward_mode not in (None, BackwardMode.NEGATE):
            return self._similarity(self._v[neigh_ixs, :], -1 * W, self._softmax_scale)

        v = self._v[ix]
        if np.all(v == 0):
            # TODO: make a method
            return _get_probs_for_zero_vec(len(neigh_ixs))

        if self._backward_mode == BackwardMode.NEGATE:
            v *= -1.0

        return self._similarity(v[None, :], W, self._softmax_scale)


class Stochastic(TmatCalculator):
    def __init__(
        self,
        conn: spmatrix,
        x: np.ndarray,
        exp: np.ndarray,
        var: np.ndarray,
        **kwargs: Any,
    ):
        super().__init__(conn, x, exp, **kwargs)
        if not hasattr(self._similarity, "hessian"):
            raise AttributeError("Similarity scheme doesn't have a `hessian` function.")
        self._var = var

    def _compute_transtion_matrix(
        self, ix: int, nbhs_ixs: np.ndarray
    ) -> Union[np.ndarray, spmatrix]:
        v = self._v[ix]
        n_neigh, n_feat = len(nbhs_ixs), self._x.shape[1]

        if np.all(v == 0):
            return _get_probs_for_zero_vec(n_neigh)
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

        # combine both to give the second order Taylor approximation. Can sometimes be negative because we
        # neglected higher order terms, so force it to be non-negative
        p = np.clip(p_0 + p_2, a_min=0, a_max=1)

        mask = np.isnan(p)
        p[mask] = 0
        c[mask] = 0

        if np.all(p == 0):
            return _get_probs_for_zero_vec(n_neigh)

        sum_ = fsum(p)
        if not np.isclose(sum_, 1.0, rtol=_RTOL):
            p[~mask] = p[~mask] / sum_

        return p, c


class Markov(TmatCalculator):
    def __init__(
        self,
        conn: spmatrix,
        x: np.ndarray,
        exp: np.ndarray,
        var: np.ndarray,
        n_samples: int = 1,
        **kwargs: Any,
    ):
        super().__init__(conn, x, exp, **kwargs)
        self._var = var
        self._n_samples = n_samples

    def _compute_transtion_matrix(
        self, ix: int, nbhs_ixs: np.ndarray
    ) -> Union[np.ndarray, spmatrix]:
        n_neigh = len(nbhs_ixs)
        W = self._x[nbhs_ixs, :] - self._x[ix, :]

        # much faster (1.8x than sampling only 1 if average)
        samples = _random_normal(
            self._x[ix],
            self._var[ix],
            n_samples=self._n_samples,
        )

        probs = np.empty((self._n_samples, n_neigh), dtype=np.float64)
        cors = np.empty((self._n_samples, n_neigh), dtype=np.float64)
        for j in prange(self._n_samples):
            probs[j], cors[j] = self._similarity(
                np.atleast_2d(samples[j]), W, self._softmax_scale
            )

        return np_mean(probs, axis=0), np_mean(cors, axis=0)
