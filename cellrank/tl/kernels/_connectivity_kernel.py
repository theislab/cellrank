from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl.kernels._bk import UnidirectionalKernel
from cellrank.tl.kernels._mixins import ConnectivityMixin

__all__ = ("ConnectivityKernel",)


@d.dedent
class ConnectivityKernel(ConnectivityMixin, UnidirectionalKernel):
    """
    Kernel which computes transition probabilities based on similarities among cells.

    As a measure of similarity, we currently support:

        - transcriptomic similarities, computed using e.g. :func:`scanpy.pp.neighbors`, see :cite:`wolf:18`.
        - spatial similarities, computed using e.g. :func:`squidpy.gr.spatial_neighbors`, see :cite:`palla:21`.

    The resulting transition matrix is symmetric and thus cannot be used to learn about the direction of the biological
    process. To include this direction, consider combining with a velocity-derived transition matrix via
    :class:`cellrank.kernels.VelocityKernel`.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    conn_key
        Key in :attr:`anndata.AnnData.obsp` where connectivity matrix describing cell-cell similarity is stored.
    check_connectivity
        Check whether the underlying KNN graph is connected.
    """

    # TODO(michalk8): chk conn docrep
    def __init__(
        self,
        adata: AnnData,
        conn_key: str = "connectivities",
        check_connectivity: bool = False,
    ):
        super().__init__(
            adata,
            conn_key=conn_key,
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
        Self and updates :attr:`transition_matrix` and :attr:`params`.
        """
        # fmt: off
        start = logg.info(f"Computing transition matrix based on `adata.obsp[{self._conn_key!r}]`")
        if self._reuse_cache({"dnorm": density_normalize, "key": self._conn_key}, time=start):
            return self

        self.transition_matrix = self._density_normalize(self._conn) if density_normalize else self._conn
        logg.info("    Finish", time=start)
        # fmt: on

        return self
