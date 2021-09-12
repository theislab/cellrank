from copy import copy

from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl.kernels import Kernel


@d.dedent
class ConnectivityKernel(Kernel):
    """
    Kernel which computes transition probabilities based on similarities among cells.

    As a measure of similarity, we currently support:

        - transcriptomic similarities, computed using e.g. :func:`scanpy.pp.neighbors`, see :cite:`wolf:18`.
        - spatial similarities, computed using e.g. :func:`squidpy.gr.spatial_neighbors`, see :cite:`palla:21`.

    The resulting transition matrix is symmetric and thus cannot be used to learn about the direction of the biological
    process. To include this direction, consider combining with a velocity-derived transition matrix via
    :class:`cellrank.tl.kernels.VelocityKernel`.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    conn_key
        Key in :attr:`anndata.AnnData.obsp` to obtain the connectivity matrix describing cell-cell similarity.
    %(cond_num)s
    check_connectivity
        Check whether the underlying KNN graph is connected.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        conn_key: str = "connectivities",
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
            conn_key=conn_key,
        )
        self._key = conn_key

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
        Self and updated :attr:`transition_matrix`.
        """

        # fmt: off
        start = logg.info(f"Computing transition matrix based on `adata.obsp[{self._key!r}]`")

        if self._reuse_cache({"dnorm": density_normalize, "key": self._key}, time=start):
            return self

        self._compute_transition_matrix(matrix=self._conn.copy(), density_normalize=density_normalize)
        logg.info("    Finish", time=start)
        # fmt: on

        return self

    def copy(self) -> "ConnectivityKernel":
        """Return a copy of self."""
        ck = ConnectivityKernel(self.adata, backward=self.backward)
        ck._params = copy(self.params)
        ck._cond_num = self.condition_number
        ck._transition_matrix = copy(self._transition_matrix)

        return ck
