# -*- coding: utf-8 -*-
"""Velocity kernel module."""
from copy import copy
from typing import Iterable, Optional
from functools import partial

from scvelo.preprocessing.moments import get_moments

import numpy as np
from cellrank import logging as logg
from scipy.sparse import issparse, csr_matrix
from cellrank.utils._docs import d
from cellrank.tools._utils import _vals_to_csr, _predict_transition_probabilities
from cellrank.tools.kernels import Kernel
from cellrank.tools._constants import Direction
from cellrank.tools.kernels._base_kernel import (
    _LOG_USING_CACHE,
    _ERROR_EMPTY_CACHE_MSG,
    AnnData,
)


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

    def compute_transition_matrix(
        self,
        backward_mode: str = "transpose",
        softmax_scale: float = 4.0,
        mode: str = "deterministic",
        seed: Optional[int] = None,
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
        seed
            Set the seed for random state, only relevant for `mode='sampling'`.

        Returns
        -------
        None
            Nothing, just makes available the following fields:

                - :paramref:`transition_matrix`
                - :paramref:`pearson_correlations`
        """

        if mode not in ["stochastic", "deterministic", "sampling"]:
            raise NotImplementedError(f"Mode `{mode}` is not implemented.")

        if mode != "deterministic" and self._direction == Direction.BACKWARD:
            logg.warning(
                f"Mode `{mode}` is currently not supported for the backward process. "
                f"Defaulting to mode `'deterministic'`"
            )
            mode = "deterministic"

        start = logg.info("Computing transition matrix based on velocity correlations")

        params = dict(  # noqa
            bwd_mode=backward_mode if self._direction == Direction.BACKWARD else None,
            sigma_corr=softmax_scale,
            mode=mode,
            random_state=seed,
        )

        # check whether we already computed such a transition matrix. If yes, load from cache
        if params == self._params:
            assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE, time=start)
            logg.info("    Finish", time=start)
            return self

        self._params = params

        # compute first and second order moments to model the distribution of the velocity vector
        if mode in ["stochastic", "sampling"]:
            np.random.seed(seed)
            velocity_expectation = get_moments(
                self.adata, self._velocity, second_order=False
            )
            velocity_variance = get_moments(
                self.adata, self._velocity, second_order=True
            )

        # define a function to compute hessian matrices
        if mode == "stochastic":
            try:
                from jax import jit, jacfwd, jacrev

                get_hessian_fwd = jit(
                    jacfwd(
                        jacrev(
                            partial(
                                _predict_transition_probabilities,
                                return_pearson_correlation=False,
                                use_jax=True,
                            )
                        )
                    )
                )
            except ImportError:
                logg.warning("TODO: install jax")
                mode = "sampling"

        # sort cells by their number of neighbors - this makes jitting more efficient
        n_neighbors = np.array((self._conn != 0).sum(1)).flatten()
        cell_ixs = np.argsort(n_neighbors)[::-1]

        # loop over all cells
        probs_list, corrs_list, rows, cols, n_obs = (
            [],
            [],
            [],
            [],
            self._gene_expression.shape[0],
        )
        for cell_ix in cell_ixs:

            # get the neighbors
            nbhs_ixs = self._conn[cell_ix, :].indices

            # get the displacement matrix. Changing dimensions b/c varying numbers of neighbors slow down autograd
            W = self._gene_expression[nbhs_ixs, :] - self._gene_expression[cell_ix, :]

            if mode == "deterministic":
                # treat `v_i` as deterministic, with no error

                if self._direction == Direction.FORWARD or (
                    self._direction == Direction.BACKWARD and backward_mode == "negate"
                ):
                    # evaluate the prediction at the actual velocity vector
                    v_i = self._velocity[cell_ix, :]

                    # for the transpose backward mode, just flip the velocity vector
                    if self._direction == Direction.BACKWARD:
                        v_i *= -1

                    if (v_i == 0).all():
                        logg.warning(f"Cell `{cell_ix}` has a zero velocity vector")
                        probs, corrs = (
                            np.ones(len(nbhs_ixs)) / len(nbhs_ixs),
                            np.zeros(len(nbhs_ixs)),
                        )
                    else:
                        probs, corrs = _predict_transition_probabilities(
                            v_i, W, softmax_scale=softmax_scale,
                        )
                else:
                    # compute how likely all neighbors are to transition to this cell
                    V = self._velocity[nbhs_ixs, :]
                    probs, corrs = _predict_transition_probabilities(
                        V, -1 * W, softmax_scale=softmax_scale,
                    )

            elif mode == "stochastic":
                # treat `v_i` as random variable and use analytical approximation to propagate the distribution

                # get the expected velocity vector
                v_i = velocity_expectation[cell_ix, :]

                if (v_i == 0).all():
                    logg.warning(f"Cell `{cell_ix}` has a zero velocity vector")
                    probs, corrs = (
                        np.ones(len(nbhs_ixs)) / len(nbhs_ixs),
                        np.zeros(len(nbhs_ixs)),
                    )
                else:
                    # get the variance of the distribution over velocity vectors
                    variances = velocity_variance[cell_ix, :]

                    # compute the Hessian tensor, and turn it into a matrix that has the diagonal elements in its rows
                    H = get_hessian_fwd(v_i, W, softmax_scale)
                    H_diag = np.concatenate([np.diag(h)[None, :] for h in H])

                    # compute zero order term
                    p_0, corrs = _predict_transition_probabilities(
                        v_i, W, softmax_scale=softmax_scale,
                    )

                    # compute second order term (note that the first order term cancels)
                    p_2 = 0.5 * H_diag.dot(variances)

                    # combine both to give the second order Taylor approximation. Can sometimes be negative because we
                    # neglected higher order terms, so force it to be non-negative
                    probs = np.where((p_0 + p_2) >= 0, p_0 + p_2, 0)

            elif mode == "sampling":
                # treat `v_i` as random variable and use Monte Carlo approximation to propagate the distribution

                # sample a single velocity vector from the distribution
                v_i = np.random.normal(
                    velocity_expectation[cell_ix, :], velocity_variance[cell_ix, :]
                )

                # use this sample to compute transition probabilities
                probs, corrs = _predict_transition_probabilities(
                    v_i, W, softmax_scale=softmax_scale,
                )

            # add the computed transition probabilities and correlations to lists
            corrs_list.extend(corrs)
            probs_list.extend(probs)
            rows.extend(np.ones(len(nbhs_ixs)) * cell_ix)
            cols.extend(nbhs_ixs)

        # collect everything in a matrix of correlations and one of transition probabilities
        probs_matrix, corrs_matrix = np.hstack(probs_list), np.hstack(corrs_list)
        probs_matrix[np.isnan(probs_matrix)] = 0
        corrs_matrix[np.isnan(corrs_matrix)] = 0

        # transform to sparse matrices
        corrs_matrix = _vals_to_csr(corrs_matrix, rows, cols, shape=(n_obs, n_obs))
        probs_matrix = _vals_to_csr(probs_matrix, rows, cols, shape=(n_obs, n_obs))

        self._compute_transition_matrix(probs_matrix, density_normalize=False)
        self._pearson_correlations = corrs_matrix

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
