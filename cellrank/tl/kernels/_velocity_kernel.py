# -*- coding: utf-8 -*-
"""Velocity kernel module."""
import os
from sys import version_info
from copy import copy, deepcopy
from math import fsum
from shutil import rmtree
from typing import Any, Union, Callable, Iterable, Optional
from tempfile import mkdtemp

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
    _HAS_JAX,
    prange,
    np_mean,
    _filter_kwargs,
    _random_normal,
    _calculate_starts,
    _reconstruct_matrices,
    _get_probs_for_zero_vec,
    _predict_transition_probabilities_jax_H,
    _predict_transition_probabilities_numpy,
)
from cellrank.tl.kernels._base_kernel import (
    TOL,
    _LOG_USING_CACHE,
    _ERROR_EMPTY_CACHE_MSG,
    AnnData,
)


class VelocityMode(ModeEnum):  # noqa
    DETERMINISTIC = "deterministic"
    STOCHASTIC = "stochastic"
    SAMPLING = "sampling"
    MONTE_CARLO = "monte_carlo"
    PROPAGATION = "propagation"


class BackwardMode(ModeEnum):  # noqa
    TRANSPOSE = "transpose"
    NEGATE = "negate"


MAX_N_ELEMS = 1024 * 1024 * 1024  # 1 million, or 8GiB in float64
TMP_DIR = None


@d.dedent
class VelocityKernel(Kernel):
    """
    Kernel which computes a transition matrix based on velocity correlations.

    This borrows ideas from both [Manno18]_ and [Bergen19]_. In short, for each cell *i*, we compute transition
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
        self._pearson_correlations = None

        self._pcors = None
        self._tmp_dir = None

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

        # filter both the velocities and the gene expression profiles to the gene subset, densify the matrices
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

    @inject_docs(m=VelocityMode, b=BackwardMode)  # don't swap the order
    @d.dedent
    def compute_transition_matrix(
        self,
        mode: str = VelocityMode.DETERMINISTIC.s,
        backward_mode: str = BackwardMode.TRANSPOSE.s,
        softmax_scale: Optional[float] = None,
        n_samples: int = 1000,
        lazy: Optional[bool] = None,
        seed: Optional[int] = None,
        **kwargs,
    ) -> "VelocityKernel":
        """
        Compute transition matrix based on velocity directions on the local manifold.

        For each cell, infer transition probabilities based on the correlation of the cell's
        velocity-extrapolated cell state with cell states of its *K* nearest neighbors.

        Parameters
        ----------
        %(velocity_mode)s
        %(velocity_backward_mode)s
        %(softmax_scale)s
        n_samples
            Number of bootstrap samples when ``mode={m.MONTE_CARLO.s!r}`` or ``mode={m.PROPAGATION.s!r}``.
        lazy
            Whether to memory-map the transition matrices when ``mode={m.PROPAGATION.s!r}``.
            If `None`, the matrices will be memory-mapped if they will be larger in total than 8 GiB.
        seed
            Set the seed for random state when the method requires ``n_samples``.
        %(parallel)s

        Returns
        -------
        :class:`cellrank.tl.kernels.VelocityKernel`
            Makes available the following fields:

                - :paramref:`transition_matrix`.
                - :paramref:`pearson_correlations`.

            If ``mode={m.PROPAGATION.s!r}``, makes also available:

                - :paramref:`transition_matrices` - tuple of length ``n_samples`` of transition matrices.
                - :paramref:`_pcors` - tuple of length ``n_samples`` of pearson correlations.
        """

        mode = VelocityMode(mode)
        backward_mode = BackwardMode(backward_mode)

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
            f"Computing transition matrix based on velocity correlations using `{mode.s!r}` mode"
        )

        if seed is None:
            seed = np.random.randint(0, 2 ** 16)
        params = dict(softmax_scale=softmax_scale, mode=mode, seed=seed)  # noqa
        if self.backward:
            params["bwd_mode"] = backward_mode.s

        # check whether we already computed such a transition matrix. If yes, load from cache
        if params == self._params:
            assert self._transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE, time=start)
            logg.info("    Finish", time=start)
            return self

        self._params = params
        self._current_ix = 0

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

        if self._tmp_dir is not None:
            rmtree(self._tmp_dir, ignore_errors=True)
            self._tmp_dir = None

        if mode == VelocityMode.PROPAGATION:
            n_elems = self._conn.nnz
            if lazy is None:
                lazy = n_samples * n_elems >= MAX_N_ELEMS
            if lazy:
                self._tmp_dir = mkdtemp(prefix="cellrank_vk_", dir=TMP_DIR)
                logg.debug(
                    f"Memory mapping sampled transition matrices and Pearson correlations under `{self._tmp_dir}`"
                )

                self._tmats = np.memmap(
                    os.path.join(self._tmp_dir, "transition_matrices.dat"),
                    mode="w+",
                    dtype=np.float64,
                    shape=(n_samples, n_elems),
                )
                self._pcors = np.memmap(
                    os.path.join(self._tmp_dir, "pearson_correlations.dat"),
                    mode="w+",
                    dtype=np.float64,
                    shape=(n_samples, n_elems),
                )
                kwargs["probs"] = self._tmats
                kwargs["cors"] = self._pcors

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
            seed=seed,
            backend=backend,
            **kwargs,
        )

        # this is ok, since the memory-mapped version returns only 1 matrix
        if isinstance(tmat, (tuple, list)):
            self._tmats, self._pcors = tuple(tmat), tuple(cmat)
            tmat, cmat = tmat[self._current_ix], cmat[self._current_ix]
        elif isinstance(self._tmats, np.memmap) and isinstance(self._pcors, np.memmap):
            logg.debug("Flushing memory mapped files")
            self._tmats.flush()
            self._pcors.flush()

            self._tmats = np.memmap(
                self._tmats.filename,
                dtype=self._tmats.dtype,
                mode="r",
                shape=self._tmats.shape,
            )
            self._pcors = np.memmap(
                self._pcors.filename,
                dtype=self._pcors.dtype,
                mode="r",
                shape=self._pcors.shape,
            )

        self._compute_transition_matrix(tmat, density_normalize=False)
        self._pearson_correlations = cmat

        logg.info("    Finish", time=start)

        return self

    @property
    def pearson_correlations(self) -> csr_matrix:  # noqa
        """The matrix of Pearson correlations."""
        return self._pearson_correlations

    @d.dedent
    def copy(self) -> "VelocityKernel":  # noqa
        """%(copy)s"""
        if isinstance(self.transition_matrices, np.memmap) or isinstance(
            self._pcors, np.memmap
        ):
            raise TypeError(
                "Transition matrices or Pearson correlations are memory mapped, copying will move them to memory."
            )

        vk = VelocityKernel(
            self.adata,
            backward=self.backward,
            vkey=self._vkey,
            xkey=self._xkey,
            gene_subset=self._gene_subset,
        )
        vk._params = copy(self.params)
        vk._cond_num = self.condition_number
        self._copy_transition_matrix(vk)
        vk._pearson_correlations = copy(self.pearson_correlations)
        vk._pcors = deepcopy(self._pcors)

        return vk

    def __del__(self):
        if self._tmp_dir is not None:
            rmtree(self._tmp_dir, ignore_errors=True)


@valuedispatch
def _dispatch_computation(mode, *_args, **_kwargs):
    raise NotImplementedError(mode)


def _run_in_parallel(fn: Callable, conn: csr_matrix, **kwargs) -> Any:
    def extractor(res):
        if res[0] is None:
            probs, cors = kwargs.get("probs", None), kwargs.get("cors", None)
            assert isinstance(
                probs, np.memmap
            ), f"Expected memory mapped data, got `{type(probs).__name__}`."
            assert isinstance(
                cors, np.memmap
            ), f"Expected memory mapped data, got `{type(cors).__name__}`."
            # get the first sample
            res = np.stack([probs[0], cors[0]])
        else:
            assert res[0].ndim in (2, 3), f"Dimension mismatch: `{res[0].ndim}`."
            res = np.concatenate(res, axis=-1)
            if res.shape[0] == 1:
                # sampling, MC, propagation with only 1
                res = res[0]

        return _reconstruct_matrices(res, conn, ixs, n_jobs=kwargs.pop("n_jobs", 1))

    fname = fn.__name__
    if fname == "_run_stochastic":
        if not _HAS_JAX:
            raise RuntimeError(
                "Install `jax` and `jaxlib` as `pip install jax jaxlib`."
            )
        ixs = np.argsort(np.array((conn != 0).sum(1)).ravel())[
            ::-1
        ]  # important - keep the descending order
    else:
        # very important for 2 reasons:
        # 1. the function requires sorted indices for correct offset calculation (we could pass it, but this is elegant)
        # 2. we get more cache hits, since we're not writing randomly into the array
        # the below case doesn't matter, because there we create an array for each process separately whereas in the
        # memory mapped case, we share a pointer to a possibly large file
        # the downside is, of course, possible uneven split between the data
        ixs = np.arange(conn.shape[0])
        if kwargs.get("probs", None) is None:
            np.random.shuffle(ixs)

    kwargs["indices"] = conn.indices
    kwargs["indptr"] = conn.indptr
    kwargs["average"] = kwargs.get("average", True)

    unit = (
        "sample" if (fname == "_run_mc") and kwargs.get("n_samples", 1) > 1 else "cell"
    )

    return parallelize(
        fn,
        ixs,
        as_array=False,
        extractor=extractor,
        unit=unit,
        **_filter_kwargs(parallelize, **kwargs),
    )(**_filter_kwargs(fn, **kwargs))


def _run_deterministic(
    ixs: Union[prange, np.ndarray],
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
        assert n_neigh, "Cell does not have any neighbors."

        W = expression[nbhs_ixs, :] - expression[ix, :]

        if not backward or (backward and backward_mode == BackwardMode.NEGATE):
            v = np.expand_dims(velocity[ix], 0)

            # for the negate backward mode
            if backward:
                v *= -1

            if np.all(v == 0):
                p, c = _get_probs_for_zero_vec(len(nbhs_ixs))
            else:
                p, c = _predict_transition_probabilities_numpy(
                    v,
                    W,
                    softmax_scale=softmax_scale,
                )
        else:
            # compute how likely all neighbors are to transition to this cell
            V = velocity[nbhs_ixs, :]
            p, c = _predict_transition_probabilities_numpy(
                V,
                -1 * W,
                softmax_scale=softmax_scale,
            )

        probs_cors[0, starts[i] : starts[i] + n_neigh] = p
        probs_cors[1, starts[i] : starts[i] + n_neigh] = c

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return probs_cors


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
    probs: Optional[np.memmap] = None,
    cors: Optional[np.memmap] = None,
) -> Optional[np.ndarray]:

    starts = _calculate_starts(indptr, ixs)
    if probs is not None and cors is not None:
        assert not average, "Memory mapping is only supported when not averaging."
        probs_cors = [probs, cors]
        offset = indptr[np.min(ixs)]  # where we start
    else:
        probs_cors = (
            np.empty((2, starts[-1]))
            if average
            else np.empty((2, n_samples, starts[-1]))
        )
        offset = 0

    for i in prange(len(ixs)):
        ix = ixs[i]
        start, end = indptr[ix], indptr[ix + 1]

        np.random.seed(seed + ix)

        nbhs_ixs = indices[start:end]
        n_neigh = len(nbhs_ixs)
        assert n_neigh, f"Cell with index `{ix}` does not have any neighbors."

        # get the displacement matrix. Changing dimensions b/c varying numbers of neighbors slow down autograd
        W = expression[nbhs_ixs, :] - expression[ix, :]

        # much faster (1.8x than sampling only 1 if average)
        samples = _random_normal(
            expectation[ix],
            variance[ix],
            n_samples=n_samples,
        )
        if average:
            probs_cors_tmp = np.empty((2, n_samples, n_neigh))
            for j in range(n_samples):
                prob, cor = _predict_transition_probabilities_numpy(
                    np.atleast_2d(samples[j]), W, softmax_scale=softmax_scale
                )
                probs_cors_tmp[0, j, :] = prob
                probs_cors_tmp[1, j, :] = cor

            probs_cors[0, starts[i] : starts[i] + n_neigh] = np_mean(
                probs_cors_tmp[0], axis=0
            )
            probs_cors[1, starts[i] : starts[i] + n_neigh] = np_mean(
                probs_cors_tmp[1], axis=0
            )
        else:
            for j in range(samples.shape[0]):
                prob, cor = _predict_transition_probabilities_numpy(
                    np.atleast_2d(samples[j]), W, softmax_scale=softmax_scale
                )
                probs_cors[0][
                    j, offset + starts[i] : offset + starts[i] + n_neigh
                ] = prob
                probs_cors[1][
                    j, offset + starts[i] : offset + starts[i] + n_neigh
                ] = cor

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    if probs is not None and cors is not None:
        return None

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
    starts = _calculate_starts(indptr, ixs)
    probs_cors = np.empty((2, starts[-1]))

    max_n = np.max(np.diff(starts))

    for i, ix in enumerate(ixs):
        start, end = indptr[ix], indptr[ix + 1]

        nbhs_ixs = indices[start:end]
        n_neigh = len(nbhs_ixs)
        assert n_neigh, "Cell does not have any neighbors."

        v = expectation[ix]

        if np.all(v == 0):
            p, c = _get_probs_for_zero_vec(len(nbhs_ixs))
        else:
            # compute the Hessian tensor, and turn it into a matrix that has the diagonal elements in its rows
            W = expression[nbhs_ixs, :] - expression[ix, :]

            W_size = W.shape[0]
            if max_n - W_size >= 10:
                max_n -= 10
            W_padded = np.pad(W, [(0, max_n - W.shape[0]), (0, 0)])

            H = _predict_transition_probabilities_jax_H(
                v, W_padded, softmax_scale=softmax_scale
            )[:W_size]
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

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return probs_cors


@_dispatch_computation.register(VelocityMode.SAMPLING)
def _(**kwargs):
    # cant 't use partial - kwargs has the attributes (could filter in principle)
    kwargs["n_samples"] = 1
    kwargs["average"] = False
    return _run_in_parallel(_run_mc, **kwargs)


@_dispatch_computation.register(VelocityMode.PROPAGATION)
def _(**kwargs):
    kwargs["average"] = False
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
