# -*- coding: utf-8 -*-
"""Velocity kernel module."""
from sys import version_info
from copy import copy
from math import fsum
from typing import Any, Union, Callable, Iterable, Optional

from scvelo.preprocessing.moments import get_moments

import numpy as np
from scipy.sparse import issparse, csr_matrix

from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.ul._utils import valuedispatch
from cellrank.tl.kernels import Kernel
from cellrank.tl._constants import _DEFAULT_BACKEND, ModeEnum
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
from cellrank.tl.kernels._base_kernel import (
    _RTOL,
    _LOG_USING_CACHE,
    _ERROR_EMPTY_CACHE_MSG,
    AnnData,
)
from cellrank.tl.kernels._velocity_schemes import Scheme, _get_scheme


class VelocityMode(ModeEnum):  # noqa
    DETERMINISTIC = "deterministic"
    STOCHASTIC = "stochastic"
    SAMPLING = "sampling"
    MONTE_CARLO = "monte_carlo"


class BackwardMode(ModeEnum):  # noqa
    TRANSPOSE = "transpose"
    NEGATE = "negate"


@d.dedent
class VelocityKernel(Kernel):
    """
    Kernel which computes a transition matrix based on RNA velocity.

    This borrows ideas from both [Manno18]_ and [Bergen20]_. In short, for each cell *i*, we compute transition
    probabilities :math:`p_{i, j}` to each cell *j* in the neighborhood of *i*. The transition probabilities are
    computed as a multinomial logistic regression where the weights :math:`w_j` (for all *j*) are given by the vector
    that connects cell *i* with cell *j* in gene expression space, and the features :math:`x_i` are given
    by the velocity vector :math:`v_i` of cell *i*.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    vkey
        Key in :paramref:`adata` ``.uns`` where the velocities are stored.
    xkey
        Key in :paramref:`adata` ``.layers`` where expected gene expression counts are stored.
    gene_subset
        List of genes to be used to compute transition probabilities.
        By default, genes from :paramref:`adata` ``.var['velocity_genes']`` are used.
    compute_cond_num
        Whether to compute condition number of the transition matrix. Note that this might be costly,
        since it does not use sparse implementation.
    check_connectivity
        Check whether the underlying KNN graph is connected.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        vkey: str = "velocity",
        xkey: str = "Ms",
        gene_subset: Optional[Iterable] = None,
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            vkey=vkey,
            xkey=xkey,
            gene_subset=gene_subset,
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
        )
        self._vkey = vkey  # for copy
        self._xkey = xkey
        self._gene_subset = gene_subset
        self._logits = None

    def _read_from_adata(self, **kwargs):
        super()._read_from_adata(**kwargs)

        # check whether velocities have been computed
        vkey = kwargs.pop("vkey", "velocity")
        if vkey not in self.adata.layers.keys():
            raise KeyError("Compute RNA velocity first as `scv.tl.velocity()`.")

        # restrict genes to a subset, i.e. velocity genes or user provided list
        gene_subset = kwargs.pop("gene_subset", None)
        subset = np.ones(self.adata.n_vars, bool)
        if gene_subset is not None:
            var_names_subset = self.adata.var_names.isin(gene_subset)
            subset &= var_names_subset if len(var_names_subset) > 0 else gene_subset
        elif f"{vkey}_genes" in self.adata.var.keys():
            subset &= np.array(self.adata.var[f"{vkey}_genes"].values, dtype=bool)

        # chose data representation to use for transcriptomic displacements
        xkey = kwargs.pop("xkey", "Ms")
        xkey = xkey if xkey in self.adata.layers.keys() else "spliced"

        # filter both the velocities and the gene expression profiles to the gene subset. Densify the matrices.
        X = np.array(
            self.adata.layers[xkey].A[:, subset]
            if issparse(self.adata.layers[xkey])
            else self.adata.layers[xkey][:, subset]
        )
        V = np.array(
            self.adata.layers[vkey].A[:, subset]
            if issparse(self.adata.layers[vkey])
            else self.adata.layers[vkey][:, subset]
        )

        # remove genes that have any Nan values (in both X and V)
        nans = np.isnan(np.sum(V, axis=0))
        if np.any(nans):
            X = X[:, ~nans]
            V = V[:, ~nans]

        # check the velocity parameters
        par_key = f"{vkey}_params"
        if par_key in self.adata.uns.keys():
            velocity_params = self.adata.uns[par_key]
        else:
            velocity_params = None
            logg.debug(
                f"Unable to load velocity parameters from `adata.uns[{par_key!r}]`"
            )

        # add to self
        self._velocity = V.astype(np.float64)
        self._gene_expression = X.astype(np.float64)
        self._velocity_params = velocity_params

    @inject_docs(m=VelocityMode, b=BackwardMode, s=Scheme)  # don't swap the order
    @d.dedent
    def compute_transition_matrix(
        self,
        mode: str = VelocityMode.DETERMINISTIC.s,
        backward_mode: str = BackwardMode.TRANSPOSE.s,
        scheme: Union[str, Callable] = Scheme.CORRELATION.s,
        softmax_scale: Optional[float] = None,
        n_samples: int = 1000,
        seed: Optional[int] = None,
        **kwargs,
    ) -> "VelocityKernel":
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
            Number of bootstrap samples when ``mode={m.MONTE_CARLO.s!r}``.
        seed
            Set the seed for random state when the method requires ``n_samples``.
        %(parallel)s

        Returns
        -------
        :class:`cellrank.tl.kernels.VelocityKernel`
            Makes available the following fields:

                - :paramref:`transition_matrix`.
                - :paramref:`logits`.
        """
        mode = VelocityMode(mode)
        backward_mode = BackwardMode(backward_mode)

        if isinstance(scheme, str):
            scheme = _get_scheme(Scheme(scheme))

        if not callable(scheme):
            raise TypeError(
                f"Expected `scheme` to be a function, found `{type(scheme)!r}`."
            )
        elif isinstance(scheme, type):
            raise TypeError(
                f"Expected `scheme` to be a function object, found function type `{type(scheme)!r}`. "
                f"Try instantiating the type."
            )

        if self.backward and mode != VelocityMode.DETERMINISTIC:
            logg.warning(
                f"Mode `{mode.s!r}` is currently not supported for the backward process. "
                f"Defaulting to mode `{VelocityMode.DETERMINISTIC.s!r}`"
            )
            mode = VelocityMode.DETERMINISTIC

        if mode == VelocityMode.STOCHASTIC and not hasattr(scheme, "hessian"):
            logg.warning(
                f"Unable to detect a method for Hessian computation. If using predefined functions, consider "
                f"installing `jax` as `pip install jax jaxlib`.\n"
                f"Defaulting to mode `{VelocityMode.MONTE_CARLO.s!r}`"
            )
            mode = VelocityMode.MONTE_CARLO

        start = logg.info(
            f"Computing transition matrix based on logits using `{mode.s!r}` mode"
        )

        if seed is None:
            seed = np.random.randint(0, 2 ** 16)
        params = dict(  # noqa: C408
            softmax_scale=softmax_scale, mode=mode.s, seed=seed, scheme=str(scheme)
        )
        if self.backward:
            params["bwd_mode"] = backward_mode.s

        # check whether we already computed such a transition matrix
        if params == self._params:
            assert self._transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE, time=start)
            logg.info("    Finish", time=start)
            return self

        self._params = params

        # compute first and second order moments to model the distribution of the velocity vector
        np.random.seed(seed)
        velocity_expectation = get_moments(
            self.adata, self._velocity, second_order=False
        ).astype(np.float64)
        velocity_variance = get_moments(
            self.adata, self._velocity, second_order=True
        ).astype(np.float64)

        if mode == VelocityMode.MONTE_CARLO and n_samples == 1:
            logg.debug("Setting mode to sampling because `n_samples=1`")
            mode = VelocityMode.SAMPLING

        backend = kwargs.pop("backend", _DEFAULT_BACKEND)
        if version_info[:2] <= (3, 6):
            logg.warning("For Python3.6, only `'threading'` backend is supported")
            backend = "threading"
        elif mode != VelocityMode.STOCHASTIC and backend == "multiprocessing":
            # this is because on jitting and pickling (cloudpickle, used by loky, handles it correctly)
            logg.warning(
                f"Multiprocessing backend is supported only for mode `{VelocityMode.STOCHASTIC.s!r}`. "
                f"Defaulting to `{_DEFAULT_BACKEND}`"
            )
            backend = _DEFAULT_BACKEND

        if softmax_scale is None:
            logg.info(
                f"Estimating `softmax_scale` using `{VelocityMode.DETERMINISTIC.s!r}` mode"
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
            softmax_scale = 1.0 / np.median(np.abs(cmat.data))
            params["softmax_scale"] = softmax_scale
            logg.info(f"Setting `softmax_scale={softmax_scale:.4f}`")

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
        self._compute_transition_matrix(tmat, density_normalize=False)
        self._logits = cmat

        logg.info("    Finish", time=start)

        return self

    @property
    def logits(self) -> csr_matrix:  # noqa
        """Array of shape ``(n_cells, n_cells)`` containing the logits."""
        return self._logits

    @d.dedent
    def copy(self) -> "VelocityKernel":  # noqa
        """%(copy)s"""
        vk = VelocityKernel(
            self.adata,
            backward=self.backward,
            vkey=self._vkey,
            xkey=self._xkey,
            gene_subset=self._gene_subset,
        )
        vk._params = copy(self.params)
        vk._cond_num = self.condition_number
        vk._transition_matrix = copy(self._transition_matrix)
        vk._logits = copy(self.logits)

        return vk


@valuedispatch
def _dispatch_computation(mode, *_args, **_kwargs):
    raise NotImplementedError(mode)


def _run_in_parallel(fn: Callable, conn: csr_matrix, **kwargs) -> Any:
    fname = fn.__name__
    if fname == "_run_stochastic":
        ixs = np.argsort(np.array((conn != 0).sum(1)).ravel())[::-1]
    else:
        ixs = np.arange(conn.shape[0])
        np.random.shuffle(ixs)

    unit = (
        "sample" if (fname == "_run_mc") and kwargs.get("n_samples", 1) > 1 else "cell"
    )
    kwargs["indices"] = conn.indices
    kwargs["indptr"] = conn.indptr

    return parallelize(
        fn,
        ixs,
        as_array=False,
        extractor=lambda res: _reconstruct_one(np.concatenate(res, axis=-1), conn, ixs),
        unit=unit,
        **_filter_kwargs(parallelize, **kwargs),
    )(**_filter_kwargs(fn, **kwargs))


def _run_deterministic(
    ixs: Union[prange, np.ndarray],
    scheme: Callable,
    indices: np.ndarray,
    indptr: np.ndarray,
    expression: np.ndarray,
    velocity: np.ndarray,
    backward: bool = False,
    backward_mode: BackwardMode = BackwardMode.TRANSPOSE,
    softmax_scale: float = 1,
    queue=None,
) -> np.ndarray:

    starts = _calculate_starts(indptr, ixs)
    probs_cors = np.empty((2, starts[-1]))

    for i in prange(len(ixs)):
        ix = ixs[i]
        start, end = indptr[ix], indptr[ix + 1]

        nbhs_ixs = indices[start:end]
        n_neigh = len(nbhs_ixs)
        assert n_neigh, f"Cell {ix} does not have any neighbors."

        W = expression[nbhs_ixs, :] - expression[ix, :]

        if not backward or (backward and backward_mode == BackwardMode.NEGATE):
            v = np.expand_dims(velocity[ix], 0)

            # for the negate backward mode
            if backward:
                v *= -1

            if np.all(v == 0):
                p, c = _get_probs_for_zero_vec(len(nbhs_ixs))
            else:
                p, c = scheme(v, W, softmax_scale)
        else:
            # compute how likely all neighbors are to transition to this cell
            V = velocity[nbhs_ixs, :]
            p, c = scheme(V, -1 * W, softmax_scale)

        probs_cors[0, starts[i] : starts[i] + n_neigh] = p
        probs_cors[1, starts[i] : starts[i] + n_neigh] = c

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return probs_cors


def _run_mc(
    ixs: Optional[np.ndarray],
    scheme: Callable,
    indices: np.ndarray,
    indptr: np.ndarray,
    expression: np.ndarray,
    expectation: np.ndarray,
    variance: np.ndarray,
    n_samples: int = 1,
    softmax_scale: float = 1,
    seed: int = 0,
    queue=None,
) -> np.ndarray:

    starts = _calculate_starts(indptr, ixs)
    probs_cors = np.empty((2, starts[-1]))

    for i in prange(len(ixs)):
        ix = ixs[i]
        start, end = indptr[ix], indptr[ix + 1]

        np.random.seed(seed + ix)

        nbhs_ixs = indices[start:end]
        n_neigh = len(nbhs_ixs)
        assert n_neigh, f"Cell {ix} does not have any neighbors."

        # get the displacement matrix. Changing dimensions b/c varying numbers of neighbors slow down autograd
        W = expression[nbhs_ixs, :] - expression[ix, :]

        # much faster (1.8x than sampling only 1 if average)
        samples = _random_normal(
            expectation[ix],
            variance[ix],
            n_samples=n_samples,
        )

        probs_cors_tmp = np.empty((2, n_samples, n_neigh))
        for j in prange(n_samples):
            prob, cor = scheme(np.atleast_2d(samples[j]), W, softmax_scale)
            probs_cors_tmp[0, j, :] = prob
            probs_cors_tmp[1, j, :] = cor

        probs_cors[0, starts[i] : starts[i] + n_neigh] = np_mean(
            probs_cors_tmp[0], axis=0
        )
        probs_cors[1, starts[i] : starts[i] + n_neigh] = np_mean(
            probs_cors_tmp[1], axis=0
        )

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return probs_cors


def _run_stochastic(
    ixs: np.ndarray,
    scheme: Callable,
    indices: np.ndarray,
    indptr: np.ndarray,
    expression: np.ndarray,
    expectation: np.ndarray,
    variance: np.ndarray,
    softmax_scale: float = 1,
    queue=None,
) -> np.ndarray:
    if not hasattr(scheme, "hessian"):
        raise AttributeError()

    starts = _calculate_starts(indptr, ixs)
    probs_cors = np.empty((2, starts[-1]))
    n_genes = expression.shape[1]

    for i, ix in enumerate(ixs):
        start, end = indptr[ix], indptr[ix + 1]

        nbhs_ixs = indices[start:end]
        n_neigh = len(nbhs_ixs)
        assert n_neigh, f"Cell {ix} does not have any neighbors."

        v = expectation[ix]

        if np.all(v == 0):
            p, c = _get_probs_for_zero_vec(len(nbhs_ixs))
        else:
            # compute the Hessian tensor, and turn it into a matrix that has the diagonal elements in its rows
            W = expression[nbhs_ixs, :] - expression[ix, :]

            H = scheme.hessian(v, W, softmax_scale)
            if H.shape == (n_neigh, n_genes):
                H_diag = H
            elif H.shape == (n_neigh, n_genes, n_genes):
                H_diag = np.array([np.diag(h) for h in H])
            else:
                raise ValueError(
                    f"Expected full Hessian matrix of shape `{(n_neigh, n_genes, n_genes)}` "
                    f"or its diagonal of shape `{(n_neigh, n_genes)}`, found `{H.shape}`."
                )

            # compute zero order term
            p_0, c = scheme(v[None, :], W, softmax_scale)

            # compute second order term (note that the first order term cancels)
            p_2 = 0.5 * H_diag.dot(variance[ix])

            # combine both to give the second order Taylor approximation. Can sometimes be negative because we
            # neglected higher order terms, so force it to be non-negative
            p = np.clip(p_0 + p_2, a_min=0, a_max=1)

            mask = np.isnan(p)
            p[mask] = 0
            c[mask] = 0

            if np.all(p == 0):
                p, c = _get_probs_for_zero_vec(len(nbhs_ixs))

            sum_ = fsum(p)
            if not np.isclose(sum_, 1.0, rtol=_RTOL):
                p[~mask] = p[~mask] / sum_

        probs_cors[0, starts[i] : starts[i] + n_neigh] = p
        probs_cors[1, starts[i] : starts[i] + n_neigh] = c

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return probs_cors


@_dispatch_computation.register(VelocityMode.SAMPLING)
def _(**kwargs):
    # cant 't use partial - parallelize needs to know whether function has `py_func` attribute + doesn't have a __name__
    kwargs["n_samples"] = 1
    return _run_in_parallel(_run_mc, **kwargs)


_dispatch_computation.register(VelocityMode.STOCHASTIC)(
    lambda **kwargs: _run_in_parallel(_run_stochastic, **kwargs)
)
_dispatch_computation.register(VelocityMode.DETERMINISTIC)(
    lambda **kwargs: _run_in_parallel(_run_deterministic, **kwargs)
)
_dispatch_computation.register(VelocityMode.MONTE_CARLO)(
    lambda **kwargs: _run_in_parallel(_run_mc, **kwargs)
)
