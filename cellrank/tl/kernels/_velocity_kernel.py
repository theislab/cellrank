from typing import Any, Optional, Sequence

from anndata import AnnData
from cellrank.ul._docs import d
from cellrank.tl.kernels._displacement_kernel import DisplacementKernel


@d.dedent
class VelocityKernel(DisplacementKernel):
    """
    Kernel which computes a transition matrix based on RNA velocity.

    This borrows ideas from both :cite:`manno:18` and :cite:`bergen:20`. In short, for each cell *i*, we compute
    transition probabilities :math:`p_{i, j}` to each cell *j* in the neighborhood of *i*. The transition probabilities
    are computed as a multinomial logistic regression where the weights :math:`w_j` (for all *j*) are given
    by the vector that connects cell *i* with cell *j* in gene expression space, and the features :math:`x_i` are given
    by the velocity vector :math:`v_i` of cell *i*.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    xkey
        Key in :attr:`anndata.AnnData.layers` where expected gene expression counts are stored.
    vkey
        Key in :attr:`anndata.AnnData.layers` where velocities are stored.
    gene_subset
        List of genes to be used to compute transition probabilities.
        If not specified, genes from :attr:`anndata.AnnData.var` ``['{vkey}_genes']`` are used.
    kwargs
        Keyword arguments for :class:`cellrank.kernels.DisplacementKernel`.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        xkey: str = "Ms",
        vkey: str = "velocity",
        gene_subset: Optional[Sequence[str]] = None,
        **kwargs: Any,
    ):
        super().__init__(
            adata,
            backward=backward,
            xkey=xkey,
            vkey=vkey,
            subset=gene_subset,
            **kwargs,
        )
