# -*- coding: utf-8 -*-
"""Velocity kernel module."""
from copy import copy, deepcopy
from math import fsum
from typing import Any, Union, Callable, Iterable, Optional

import numpy as np
from scipy.sparse import issparse, csr_matrix

from scvelo.preprocessing.moments import get_moments

from cellrank import logging as logg
from cellrank.utils._docs import d, inject_docs
from cellrank.utils._utils import valuedispatch
from cellrank.tools.kernels import Kernel
from cellrank.tools._constants import ModeEnum
from cellrank.utils._parallelize import parallelize
from cellrank.tools.kernels._utils import (
    _HAS_JAX,
    njit,
    prange,
    jit_kwargs,
    _filter_kwargs,
    _random_normal,
    _calculate_starts,
    _reconstruct_matrices,
    _get_probs_for_zero_vec,
    _predict_transition_probabilities_jax_H,
    _predict_transition_probabilities_numpy,
)
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
    PROPAGATION = "propagation"


@d.dedent
class VelocityKernel(Kernel):
    """
    Kernel which computes a transition matrix based on velocity correlations.

    This borrows ideas from both [Manno18]_ and [Bergen19]_. In short, for each cell *i*, we compute transition
    probabilities :math:`p_{i, j}` to each cell *j* in the neighborhood of *i*. The transition probabilities are
    computed as a multinominal logistic regression where the weights :math:`w_j` (for all *j*) are given by the vector
    that connects cell *i* with cell *j* in gene expression space, and the features :math:`x_i` are given
    by the velocity vector :math:`v_i` of cell *i*.

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

        self._tmats = None
        self._pcors = None

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
                f"Unable to load velocity parameters from `adata.uns[{par_key}!]`"
            )

        # add to self
        self._velocity = V.astype(np.float64)
        self._gene_expression = X.astype(np.float64)
        self._velocity_params = velocity_params

    @d.dedent
    @inject_docs(m=VelocityMode)
    def compute_transition_matrix(
        self,
        mode: str = VelocityMode.DETERMINISTIC.s,
        backward_mode: str = "transpose",
        softmax_scale: float = 4.0,
        n_samples: int = 1000,
        seed: Optional[int] = None,
        use_numba: Optional[bool] = False,
        **kwargs,
    ) -> "VelocityKernel":
        """
        Compute transition matrix based on velocity directions on the local manifold.

        For each cell, infer transition probabilities based on the correlation of the cell's
        velocity-extrapolated cell state with cell states of its *K* nearest neighbors.

        Parameters
        ----------
        mode
            How to compute transition probabilities. Valid options are:

                - `{m.DETERMINISTIC.s!r}` - deterministic computation that doesn't propagate uncertainty
                - `{m.MONTE_CARLO.s!r}` - Monte Carlo average of randomly sampled velocity vectors
                - `{m.STOCHASTIC.s!r}` - second order approximation, only available when :mod:`jax` is installed.
                - `{m.SAMPLING.s!r}` - sample 1 transition matrix from velocity distribution
                - `{m.PROPAGATION.s!r}` - same as `{m.MONTE_CARLO.s!r}`, but does not average the vectors.
                  Instead, it saves the sampled transition matrices to :paramref:`_t_mats` to be used
                  for later uncertainty estimation. It is generally faster then `{m.MONTE_CARLO.s!r}`,
                  but also less memory efficient
        backward_mode
            Only matters if initialized as :paramref:`backward` =`True`.  Valid options are:

                - `'transpose'` - compute transitions from neighboring cells `j` to cell `i`
                - `'negate'` - negate the velocity vector
        softmax_scale
            Scaling parameter for the softmax.
        n_samples
            Number of bootstrap samples when :paramref:`mode` is `{m.MONTE_CARLO.s!r}` or `{m.PROPAGATION.s!r}`.
        seed
            Set the seed for random state when the method requires :paramref:`n_samples`.
        use_numba
            Use :mod:`numba` optimized functions. Only available if `:paramref:`mode` != `{m.STOCHASTIC.s!r}`.
            Note that this options disables the progress bar, but can be faster for.
            If `None`, it enables running :mod:`numba` jitted functions in conjunction with default parallelization
            options.
        %(parallel)s

        Returns
        -------
        :class:`cellrank.tl.kernels.VelocityKernel`
            Makes available the following fields:

                - :paramref:`transition_matrix`
                - :paramref:`pearson_correlations`

            If :paramref:`mode` `={m.PROPAGATION.s!r}`, makes also available:

                - :paramref:`_tmats` - tuple of lenth :paramref:`n_samples` of transition matrices
                - :paramref:`_pcors` - tuple of lenth :paramref:`n_samples` of Pearson correlations
        """

        mode = VelocityMode(mode)

        if self.backward and mode != VelocityMode.DETERMINISTIC:
            logg.warning(
                f"Mode `{mode.s!r}` is currently not supported for the backward process. "
                f"Defaulting to mode `{VelocityMode.DETERMINISTIC.s!r}`"
            )
            mode = VelocityMode.DETERMINISTIC
        if mode == VelocityMode.STOCHASTIC and not _HAS_JAX:
            logg.warning(
                f"Unable to detect `jax` installation. Consider installing it as `pip install jax jaxlib`.\n"
                f"Defaulting to mode `{VelocityMode.MONTE_CARLO.s!r}`"
            )
            mode = VelocityMode.MONTE_CARLO

        start = logg.info(
            f"Computing transition matrix based on velocity correlations using mode `{mode.s!r}`"
        )

        if seed is None:
            seed = np.random.randint(0, 2 ** 16)
        params = dict(  # noqa
            bwd_mode=backward_mode if self.backward else None,
            sigma_corr=softmax_scale,
            mode=mode,
            seed=seed,
        )

        # check whether we already computed such a transition matrix. If yes, load from cache
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

        tmat, cmat = _dispatch_computation(
            mode,
            conn=self._conn,
            expression=self._gene_expression,
            velocity=self._velocity,
            expectation=velocity_expectation,
            variance=velocity_variance,
            softmax_scale=softmax_scale,
            backward=self.backward,
            backward_mode=backward_mode,
            n_samples=n_samples,
            use_numba=use_numba,
            seed=seed,
            **kwargs,
        )
        if isinstance(tmat, (tuple, list)):
            self._tmats, self._pcors = tmat, cmat
            tmat, cmat = tmat[0], cmat[0]

        self._compute_transition_matrix(tmat, density_normalize=False)
        self._pearson_correlations = cmat

        logg.info("    Finish", time=start)

        return self

    @property
    def pearson_correlations(self) -> csr_matrix:  # noqa
        """The matrix of Pearson correlations."""
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
        vk._tmats = deepcopy(self._tmats)
        vk._pcors = deepcopy(self._pcors)

        return vk


@valuedispatch
def _dispatch_computation(mode, *_args, **_kwargs):
    raise NotImplementedError(mode)


def _run_in_parallel(fn: Callable, conn: csr_matrix, **kwargs) -> Any:
    def extractor(res):
        assert res[0].ndim in (2, 3), f"Dimension mismatch: `{res[0].ndim}`"

        res = np.concatenate(res, axis=-1)
        if res.shape[0] == 1:
            # sampling, MC, propagation with only 1
            res = res[0]

        return _reconstruct_matrices(
            res, conn, resort_ixs, n_jobs=kwargs.pop("n_jobs", 1)
        )

    fname = fn.__name__
    kwargs["average"] = kwargs.get("average", True)
    if not kwargs["average"]:
        assert fname == "_run_mc", "Invalid function."
        ixs = np.arange(conn.shape[0])
        resort_ixs = None
        unit = "sample"
    else:
        ixs = np.argsort(np.array((conn != 0).sum(1)).ravel())
        resort_ixs = ixs
        unit = (
            "sample"
            if (fname == "_run_mc") and kwargs.get("n_samples", 1) > 1
            else "cell"
        )

    kwargs["indices"] = conn.indices
    kwargs["indptr"] = conn.indptr

    return parallelize(
        fn,
        ixs,
        as_array=False,
        extractor=extractor,
        unit=unit,
        **_filter_kwargs(parallelize, **kwargs),
    )(**_filter_kwargs(fn, **kwargs))


def _run(fn: Callable, **kwargs) -> Any:

    use_numba = kwargs.pop("use_numba", False)
    conn = kwargs.pop("conn")

    kwargs["indices"] = conn.indices
    kwargs["indptr"] = conn.indptr

    if fn is _run_stochastic and not _HAS_JAX:
        raise RuntimeError("Install `jax` and `jaxlib` as `pip install jax jaxlib`.")

    if not use_numba:
        # can be disabled because of 2 reasons: numba is not installed or manually disabled
        # remove the outer numba loop
        if use_numba is not None and hasattr(fn, "py_func"):
            logg.debug(f"Disabling numba jitting for function `{fn.__name__!r}`")
            fn = fn.py_func
        return _run_in_parallel(fn, conn, **kwargs)

    kwargs["ixs"] = np.arange(conn.shape[0], dtype=np.int32)
    ks = _filter_kwargs(fn, **kwargs)

    return _reconstruct_matrices(fn(**ks), conn, ixs=None)


@njit(parallel=True, **jit_kwargs)
def _run_deterministic(
    ixs: Union[prange, np.ndarray],
    indices: np.ndarray,
    indptr: np.ndarray,
    expression: np.ndarray,
    velocity: np.ndarray,
    backward: bool = False,
    backward_mode: str = "negate",
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
        W = expression[nbhs_ixs, :] - expression[ix, :]

        if not backward or (backward and backward_mode == "negate"):
            v = np.expand_dims(velocity[ix], 0)

            # for the transpose backward mode, just flip the velocity vector
            if backward:
                v *= -1

            if np.all(v == 0):
                p, c = _get_probs_for_zero_vec(len(nbhs_ixs))
            else:
                p, c = _predict_transition_probabilities_numpy(
                    v, W, softmax_scale=softmax_scale,
                )
        else:
            # compute how likely all neighbors are to transition to this cell
            V = velocity[nbhs_ixs, :]
            p, c = _predict_transition_probabilities_numpy(
                V, -1 * W, softmax_scale=softmax_scale,
            )

        probs_cors[0, starts[i] : starts[i] + n_neigh] = p
        probs_cors[1, starts[i] : starts[i] + n_neigh] = c

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return probs_cors


@njit(parallel=True, **jit_kwargs)
def _run_mc(
    ixs: Optional[np.ndarray],
    indices: np.ndarray,
    indptr: np.ndarray,
    expression: np.ndarray,
    expectation: np.ndarray,
    variance: np.ndarray,
    n_samples: int = 1,
    average: bool = True,
    softmax_scale: float = 1,
    seed: int = 0,
    queue=None,
) -> np.ndarray:

    starts = _calculate_starts(indptr, ixs)
    if average:
        probs_cors = np.empty((1, 2, starts[-1]))
    else:
        probs_cors = np.empty((n_samples, 2, starts[-1]))

    for i in prange(len(ixs)):
        ix = ixs[i]
        np.random.seed(seed + ix)

        start, end = indptr[ix], indptr[ix + 1]
        nbhs_ixs = indices[start:end]
        n_neigh = len(nbhs_ixs)
        # get the displacement matrix. Changing dimensions b/c varying numbers of neighbors slow down autograd
        W = expression[nbhs_ixs, :] - expression[ix, :]

        if average:
            v = np.atleast_2d(
                _random_normal(
                    expectation[ix], variance[ix], n_samples=n_samples, average=True,
                )
            )

            p, c = _predict_transition_probabilities_numpy(
                v, W, softmax_scale=softmax_scale
            )

            probs_cors[0, 0, starts[i] : starts[i] + n_neigh] = p
            probs_cors[0, 1, starts[i] : starts[i] + n_neigh] = c
        else:
            samples = _random_normal(
                expectation[ix], variance[ix], n_samples=n_samples, average=False,
            )
            for j in prange(samples.shape[0]):
                prob, cor = _predict_transition_probabilities_numpy(
                    np.atleast_2d(samples[j]), W, softmax_scale=softmax_scale
                )
                probs_cors[j, 0, starts[i] : starts[i] + n_neigh] = prob
                probs_cors[j, 1, starts[i] : starts[i] + n_neigh] = cor

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return probs_cors


def _run_stochastic(
    ixs: np.ndarray,
    indices: np.ndarray,
    indptr: np.ndarray,
    expression: np.ndarray,
    expectation: np.ndarray,
    variance: np.ndarray,
    softmax_scale: float = 1,
    queue=None,
) -> np.ndarray:
    # TODO: can't use numba yet, since we rely on JAX for hessians
    starts = _calculate_starts(indptr, ixs)
    probs_cors = np.empty((2, starts[-1]))

    for i, ix in enumerate(ixs):
        start, end = indptr[ix], indptr[ix + 1]
        nbhs_ixs = indices[start:end]
        n_neigh = len(nbhs_ixs)
        v = expectation[ix]

        if np.all(v == 0):
            p, c = _get_probs_for_zero_vec(len(nbhs_ixs))
        else:
            # compute the Hessian tensor, and turn it into a matrix that has the diagonal elements in its rows
            W = expression[nbhs_ixs, :] - expression[ix, :]

            H = _predict_transition_probabilities_jax_H(
                v, W, softmax_scale=softmax_scale
            )
            H_diag = np.array([np.diag(h) for h in H])

            # compute zero order term
            p_0, c = _predict_transition_probabilities_numpy(
                v[None, :], W, softmax_scale=softmax_scale
            )

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
            if not np.isclose(sum_, 1.0, rtol=TOL):
                p[~mask] = p[~mask] / sum_

        probs_cors[0, starts[i] : starts[i] + n_neigh] = p
        probs_cors[1, starts[i] : starts[i] + n_neigh] = c

        queue.put(1)

    queue.put(None)

    return probs_cors


@_dispatch_computation.register(VelocityMode.SAMPLING)
def _(**kwargs):
    # cant 't use partial - parallelize needs to know whether function has `py_func` attribute + doesn't have a __name__
    kwargs["n_samples"] = 1
    return _run(_run_mc, **kwargs)


@_dispatch_computation.register(VelocityMode.PROPAGATION)
def _(**kwargs):
    kwargs["average"] = False
    return _run(_run_mc, **kwargs)


_dispatch_computation.register(VelocityMode.STOCHASTIC)(
    lambda **kwargs: _run_in_parallel(_run_stochastic, **kwargs)
)
_dispatch_computation.register(VelocityMode.DETERMINISTIC)(
    lambda **kwargs: _run(_run_deterministic, **kwargs)
)
_dispatch_computation.register(VelocityMode.MONTE_CARLO)(
    lambda **kwargs: _run(_run_mc, **kwargs)
)
