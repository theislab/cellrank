# -*- coding: utf-8 -*-
"""Palantir kernel module."""
from copy import copy

import numpy as np

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import _bias_knn, _connected
from cellrank.tl.kernels import Kernel
from cellrank.tl._constants import Direction
from cellrank.tl.kernels._base_kernel import (
    _LOG_USING_CACHE,
    _ERROR_EMPTY_CACHE_MSG,
    AnnData,
    _dtype,
)


@d.dedent
class PalantirKernel(Kernel):
    """
    Kernel which computes transition probabilities in a similar way to *Palantir*, see [Setty19]_.

    *Palantir* computes a KNN graph in gene expression space and a pseudotime, which it then uses to direct the edges of
    the KNN graph, such that they are more likely to point into the direction of increasing pseudotime. To avoid
    disconnecting the graph, it does not remove all edges that point into the direction of decreasing pseudotime
    but keeps the ones that point to nodes inside a close radius. This radius is chosen according to the local density.

    The implementation presented here won't exactly reproduce the original *Palantir* algorithm (see below)
    but the results are qualitatively very similar.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :paramref:`adata` ``.obs`` where the pseudotime is stored.
    compute_cond_num
        Whether to compute condition number of the transition matrix. Note that this might be costly,
        since it does not use sparse implementation.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        time_key: str = "dpt_pseudotime",
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=time_key,
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
        )
        self._time_key = time_key

    def _read_from_adata(self, **kwargs):
        super()._read_from_adata(**kwargs)

        time_key = kwargs.pop("time_key", "dpt_pseudotime")
        if time_key not in self.adata.obs.keys():
            raise KeyError(f"Could not find time key in `adata.obs[{time_key!r}]`.")

        self._pseudotime = np.array(self.adata.obs[time_key]).astype(_dtype)

        if np.any(np.isnan(self._pseudotime)):
            raise ValueError("Encountered NaN values in pseudotime.")

        logg.debug("Clipping the pseudotime to 0-1 range")
        self._pseudotime = np.clip(self._pseudotime, 0, 1)

    def compute_transition_matrix(
        self, k: int = 3, density_normalize: bool = True
    ) -> "PalantirKernel":
        """
        Compute transition matrix based on KNN graph and pseudotemporal ordering.

        This is a re-implementation of the Palantir algorithm by [Setty19]_.
        Note that this won't exactly reproduce the original Palantir results, for three reasons:

            - Palantir computes the KNN graph in a scaled space of diffusion components.
            - Palantir uses its own pseudotime to bias the KNN graph which is not implemented here.
            - Palantir uses a slightly different mechanism to ensure the graph remains connected when removing edges
              that point into the "pseudotime past".

        If you would like to reproduce the original results, please use the original Palantir algorithm.

        Parameters
        ----------
        k
            Number of neighbors to keep for each node, regardless of pseudotime.
            This is done to ensure that the graph remains connected.
        density_normalize
            Whether or not to use the underlying KNN graph for density normalization.

        Returns
        -------
        :class:`cellrank.tl.kernels.PalantirKernel`
            Makes :paramref:`transition_matrix` available.
        """

        start = logg.info("Computing transition matrix based on Palantir-like kernel")

        # get the connectivities and number of neighbors
        n_neighbors = (
            self.adata.uns.get("neighbors", {})
            .get("params", {})
            .get("n_neighbors", None)
        )
        if n_neighbors is None:
            logg.warning(
                "Could not find 'n_neighbors' in `adata.uns['neighbors']['params']`. Using an estimate"
            )
            n_neighbors = np.min(self._conn.sum(1))

        params = dict(k=k, dnorm=density_normalize, n_neighs=n_neighbors)  # noqa
        if params == self._params:
            assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE)
            logg.info("    Finish", time=start)
            return self

        self._params = params

        # handle backward case and run biasing function
        pseudotime = (
            np.max(self.pseudotime) - self.pseudotime
            if self._direction == Direction.BACKWARD
            else self.pseudotime
        )
        biased_conn = _bias_knn(
            conn=self._conn, pseudotime=pseudotime, n_neighbors=n_neighbors, k=k
        ).astype(_dtype)
        # make sure the biased graph is still connected
        if not _connected(biased_conn):
            logg.warning("Biased KNN graph is disconnected")

        self._compute_transition_matrix(
            matrix=biased_conn, density_normalize=density_normalize
        )

        logg.info("    Finish", time=start)

        return self

    @property
    def pseudotime(self) -> np.array:
        """Pseudotemporal ordering of cells."""
        return self._pseudotime

    def copy(self) -> "PalantirKernel":
        """Return a copy of self."""
        pk = PalantirKernel(self.adata, backward=self.backward, time_key=self._time_key)
        pk._pseudotime = copy(self.pseudotime)
        pk._params = copy(self._params)
        pk._cond_num = self.condition_number
        pk._transition_matrix = copy(self._transition_matrix)

        return pk

    def __invert__(self) -> "PalantirKernel":
        super().__invert__()
        self._pseudotime = 1 - self.pseudotime
        return self
