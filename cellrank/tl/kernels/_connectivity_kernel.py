"""Connectivity kernel module."""
from copy import copy

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl.kernels import Kernel
from cellrank.tl.kernels._base_kernel import (
    _LOG_USING_CACHE,
    _ERROR_EMPTY_CACHE_MSG,
    AnnData,
)


@d.dedent
class ConnectivityKernel(Kernel):
    """
    Kernel which computes transition probabilities based on transcriptomic similarities.

    As a measure for transcriptomic similarity, we use the weighted KNN graph computed using
    :func:`scanpy.pp.neighbors`, see [Wolf18]_.
    By definition, the resulting transition matrix is symmetric and cannot be used to learn about the direction of the
    developmental process under consideration.
    However, the velocity-derived transition matrix from :class:`cellrank.tl.kernels.VelocityKernel` can be combined
    with the similarity-based transition matrix as a means of regularization.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
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
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
        )

    def compute_transition_matrix(
        self, density_normalize: bool = True
    ) -> "ConnectivityKernel":
        """
        Compute transition matrix based on transcriptomic similarity.

        Uses symmetric, weighted KNN graph to compute symmetric transition matrix. The connectivities are computed
        using :func:`scanpy.pp.neighbors`. Depending on the parameters used there, they can be UMAP connectivities or
        gaussian-kernel-based connectivities with adaptive kernel width.

        Parameters
        ----------
        density_normalize
            Whether or not to use the underlying KNN graph for density normalization.

        Returns
        -------
        :class:`cellrank.tl.kernels.ConnectivityKernel`
            Makes :paramref:`transition_matrix` available.
        """

        start = logg.info("Computing transition matrix based on connectivities")

        params = {"dnorm": density_normalize}
        if params == self.params:
            assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE)
            logg.info("    Finish", time=start)
            return self

        self._params = params
        self._compute_transition_matrix(
            matrix=self._conn.copy(), density_normalize=density_normalize
        )

        logg.info("    Finish", time=start)

        return self

    def copy(self) -> "ConnectivityKernel":
        """Return a copy of self."""
        ck = ConnectivityKernel(self.adata, backward=self.backward)
        ck._params = copy(self.params)
        ck._cond_num = self.condition_number
        ck._transition_matrix = copy(self._transition_matrix)

        return ck
