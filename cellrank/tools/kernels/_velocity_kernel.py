# -*- coding: utf-8 -*-
"""Velocity kernel module."""
from copy import copy
from math import fsum
from typing import List, Iterable, Optional
from inspect import signature
from functools import partial

import numpy as np
from scipy.sparse import issparse, spmatrix, csr_matrix

from scvelo.preprocessing.moments import get_moments

from cellrank import logging as logg
from cellrank.utils._docs import d
from cellrank.tools._utils import _predict_transition_probabilities
from cellrank.utils._utils import valuedispatch
from cellrank.tools.kernels import Kernel
from cellrank.tools._constants import ModeEnum, Direction
from cellrank.utils._parallelize import parallelize
from cellrank.tools.kernels._base_kernel import (
    _LOG_USING_CACHE,
    _ERROR_EMPTY_CACHE_MSG,
    AnnData,
)

TOL = 1e-12


class VelocityMode(ModeEnum):
    """Method to calculate the transition matrix of VelocityKernel."""

    DETERMINISTIC = "deterministic"
    STOCHASTIC = "stochastic"
    SAMPLING = "sampling"
    MONTE_CARLO = "monte_carlo"


@d.dedent
class VelocityKernel(Kernel):
    """
    Kernel which computes a transition matrix based on velocity correlations.

    This borrows ideas from both [Manno18]_ and [Bergen19]_. In short, for each cell *i*, we compute transition
    probabilities :math:`p_{i, j}` to each cell *j* in the neighborhood of *i*. The transition probabilities are
    computed as a multinominal logistic regression where the weights :math:`w_j` (for all *j*) are given by the vector
    that connects cell *i* with cell *j* in gene expression space, and the features :math:`x_i` are given
    by the velocity vector :math:`v_i` of cell *i*.

    Optionally, we propagate uncertainty in the velocity vectors forward into transition probabilities using either
    an analytical approximation or a Monte Carlo approach.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    vkey
        Key in :paramref:`adata` `.uns` where the velocities are stored.
    xkey
        Key in :paramref:`adata` `.layers` where expected gene expression counts are stored.
    gene_subset
        List of genes to be used to compute transition probabilities. By default, the `velocity_genes` of
        :paramref:`adata` `. var` are used.
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
        self._pearson_correlations = None

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
        if vkey + "_params" in self.adata.uns.keys():
            velocity_params = self.adata.uns[vkey + "_params"]
        else:
            velocity_params = None
            logg.debug("Unable to load velocity parameters")

        # add to self
        self._velocity = V
        self._gene_expression = X
        self._velocity_params = velocity_params

    @d.dedent
    def compute_transition_matrix(
        self,
        backward_mode: str = "transpose",
        softmax_scale: float = 4.0,
        mode: str = "deterministic",
        seed: Optional[int] = None,
        n_samples: int = 1000,
        **kwargs,
    ) -> "VelocityKernel":
        """
        Compute transition matrix based on velocity directions on the local manifold.

        For each cell, infer transition probabilities based on the correlation of the cell's
        velocity-extrapolated cell state with cell states of its *K* nearest neighbors.

        Parameters
        ----------
        backward_mode
            Options are `['transpose', 'negate']`. Only matters if initialized as :paramref:`backward` =`True`.
        softmax_scale
            Scaling parameter for the softmax.
        mode
            How to compute transition probabilities. Options are "stochastic" (propagate uncertainty analytically),
            "deterministic" (don't propagate uncertainty) and "sampling" (sample from velocity distribution).
        n_samples
            Number of bootstrap samples when :paramref:`mode` `='monte_carlo'`.
        seed
            Set the seed for random state, only relevant for `mode='sampling'`.
        %(parallel)s

        Returns
        -------
        None
            Nothing, just makes available the following fields:

                - :paramref:`transition_matrix`
                - :paramref:`pearson_correlations`
        """

        mode = VelocityMode(mode)

        if mode != VelocityMode.DETERMINISTIC and self.backward:
            logg.warning(
                f"VelocityMode `{mode}` is currently not supported for the backward process. "
                f"Defaulting to mode `{VelocityMode.DETERMINISTIC.value!r}`"
            )
            mode = VelocityMode.DETERMINISTIC

        start = logg.info("Computing transition matrix based on velocity correlations")

        params = dict(  # noqa
            bwd_mode=backward_mode if self._direction == Direction.BACKWARD else None,
            sigma_corr=softmax_scale,
            mode=mode,
            seed=seed,
        )

        # check whether we already computed such a transition matrix. If yes, load from cache
        if params == self._params:
            assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE, time=start)
            logg.info("    Finish", time=start)
            return self

        self._params = params

        # compute first and second order moments to model the distribution of the velocity vector
        # if mode in ["stochastic", "sampling"]:
        np.random.seed(seed)
        velocity_expectation = get_moments(
            self.adata, self._velocity, second_order=False
        )
        velocity_variance = get_moments(self.adata, self._velocity, second_order=True)

        t_mat, corrs_mat = _dispatch_computation(
            VelocityMode(mode),
            conn=self._conn,
            expression=self._gene_expression,
            velocity=self._velocity,
            expectation=velocity_expectation,
            variance=velocity_variance,
            softmax_scale=softmax_scale,
            backward=self.backward,
            backward_mode=backward_mode,
            **kwargs,
        )

        self._compute_transition_matrix(t_mat, density_normalize=False)
        self._pearson_correlations = corrs_mat

        logg.info("    Finish", time=start)

        return self

    @property
    def pearson_correlations(self) -> csr_matrix:
        """TODO."""
        return self._pearson_correlations

    def copy(self) -> "VelocityKernel":
        """Return a copy of self."""
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
        vk._pearson_correlations = copy(self.pearson_correlations)

        return vk


def _get_displacement(ix, conn, expression):
    nbhs_ixs = conn[ix, :].indices

    # get the displacement matrix. Changing dimensions b/c varying numbers of neighbors slow down autograd
    W = expression[nbhs_ixs, :] - expression[ix, :]

    return W, nbhs_ixs


@valuedispatch
def _dispatch_computation(mode, *_args, **_kwargs):
    raise NotImplementedError(mode)


def _run_in_parallel(fn, conn, sort_by_neighbors: bool = False, **kwargs):
    def _to_csr_matrix(data: List[float], mat: csr_matrix):
        return csr_matrix((data, mat.indices, mat.indptr))

    if not sort_by_neighbors:
        ixs = np.arange(conn.shape[0])
    else:
        n_neighbors = np.array((conn != 0).sum(1)).flatten()
        ixs = np.argsort(n_neighbors)[::-1]

    params = signature(parallelize).parameters
    res = parallelize(
        fn,
        ixs,
        as_array=True,
        extractor=np.vstack,
        unit="cell",
        **{k: v for k, v in kwargs.items() if k in params},
    )()
    assert res.ndim == 2 and res.shape[1] == 2, "Sanity check failed."
    probs, cors = res[:, 0], res[:, 1]

    if not sort_by_neighbors:
        return _to_csr_matrix(probs, conn), _to_csr_matrix(cors, conn)

    conn = conn[ixs, :][:, ixs]  # get the correct indices and indptrs
    aixs = np.argsort(ixs)

    probs, conn = _to_csr_matrix(probs, conn), _to_csr_matrix(cors, conn)

    return probs[aixs, :][:, aixs], conn[aixs, :][:, aixs]


@_dispatch_computation.register(VelocityMode.MONTE_CARLO)
def _(
    conn: csr_matrix,
    expression: np.ndarray,
    expectation: np.ndarray,
    variance: np.ndarray,
    n_samples: int = 1,
    softmax_scale: float = 1,
    **kwargs,
):
    def runner(ixs, queue):
        probs, cors = [], []

        for ix in ixs:
            W, nbhs_ixs = _get_displacement(ix, conn, expression)
            v = _random_normal(expectation[ix], variance[ix], n_samples=n_samples)

            p, c = _predict_transition_probabilities(v, W, softmax_scale=softmax_scale,)
            probs.extend(p)
            cors.extend(c)

            queue.put(1)

        queue.put(None)

        return np.c_[probs, cors]

    return _run_in_parallel(runner, conn, **kwargs)


@_dispatch_computation.register(VelocityMode.SAMPLING)
def _(
    conn: csr_matrix,
    expression: np.ndarray,
    expectation: np.ndarray,
    variance: np.ndarray,
    n_samples: int = 1,
    softmax_scale: float = 1,
    **kwargs,
):
    loc = locals()
    loc["n_samples"] = 1
    _ = loc.pop("kwargs")
    return _dispatch_computation(VelocityMode.MONTE_CARLO, **loc, **kwargs)


@_dispatch_computation.register(VelocityMode.DETERMINISTIC)
def _(
    conn: spmatrix,
    expression: np.ndarray,
    velocity: np.ndarray,
    backward: bool = False,
    backward_mode: str = "negate",
    softmax_scale: float = 1,
    **kwargs,
):
    def runner(ixs, queue):
        probs, cors = [], []

        for ix in ixs:
            W, nbhs_ixs = _get_displacement(ix, conn, expression)

            if not backward or (backward and backward_mode == "negate"):
                # evaluate the prediction at the actual velocity vector
                v = velocity[ix]

                # for the transpose backward mode, just flip the velocity vector
                if backward:
                    v *= -1

                if np.all(v == 0):
                    p, c = _get_probs_for_zero_vec(len(nbhs_ixs))
                else:
                    p, c = _predict_transition_probabilities(
                        v, W, softmax_scale=softmax_scale,
                    )
            else:
                # compute how likely all neighbors are to transition to this cell
                V = velocity[nbhs_ixs, :]
                p, c = _predict_transition_probabilities(
                    V, -1 * W, softmax_scale=softmax_scale,
                )

            probs.extend(p)
            cors.extend(p)

            queue.put(1)

        queue.put(None)

        return np.c_[probs, cors]

    return _run_in_parallel(runner, conn, **kwargs)


@_dispatch_computation.register(VelocityMode.STOCHASTIC)
def _(
    conn: spmatrix,
    expression: np.ndarray,
    expectation: np.ndarray,
    variance: np.ndarray,
    softmax_scale: float = 1,
    **kwargs,
):
    def runner(ixs, queue):
        from jax import jit, jacfwd, jacrev

        get_hessian_fwd = jit(jacfwd(jacrev(_predict_transition_probabilities_jax)))

        probs, cors = [], []

        for ix in ixs:
            v = expectation[ix]
            W, nbhs_ixs = _get_displacement(ix, conn, expression)

            if np.all(v == 0):
                p, c = _get_probs_for_zero_vec(len(nbhs_ixs))
                sum_ = fsum(p)
                if not np.isclose(sum_, 1.0, rtol=TOL):
                    p /= sum_
            else:
                # get the variance of the distribution over velocity vectors
                variances = variance[ix]

                # compute the Hessian tensor, and turn it into a matrix that has the diagonal elements in its rows
                H = get_hessian_fwd(v, W, softmax_scale)
                H_diag = np.concatenate([np.diag(h)[None, :] for h in H])

                # compute zero order term
                p_0, c = _predict_transition_probabilities(
                    v, W, softmax_scale=softmax_scale,
                )

                # compute second order term (note that the first order term cancels)
                p_2 = 0.5 * H_diag.dot(variances)

                # combine both to give the second order Taylor approximation. Can sometimes be negative because we
                # neglected higher order terms, so force it to be non-negative
                p = np.clip(p_0 + p_2, a_min=0, a_max=1)
                mask = np.isnan(p)
                p[mask] = 0

                if np.all(p == 0):
                    p, c = _get_probs_for_zero_vec(len(nbhs_ixs))

                sum_ = fsum(p)
                if not np.isclose(sum_, 1.0, rtol=TOL):
                    p[~mask] = p[~mask] / sum_

            probs.extend(p.astype(np.float64))  # can be float32 from JAX
            cors.extend(c.astype(np.float64))

            queue.put(1)

        queue.put(None)

        return np.c_[probs, cors]

    return _run_in_parallel(runner, conn, sort_by_neighbors=True, **kwargs)


def _get_probs_for_zero_vec(size: int):
    # float32 doesn't have enough precision
    return (
        np.ones(size, dtype=np.float64) / size,
        np.zeros(size, dtype=np.float64),
    )


def _random_normal(mean: np.ndarray, var: np.ndarray, n_samples=1):
    if mean.ndim != 1:
        raise ValueError()
    if mean.shape != var.shape:
        raise ValueError()

    if n_samples == 1:
        return np.random.normal(mean, var)

    return np.mean(np.random.normal(mean, var, size=(n_samples, mean.shape[0])), axis=0)


_predict_transition_probabilities_jax = partial(
    _predict_transition_probabilities, return_pearson_correlation=False, use_jax=True,
)
