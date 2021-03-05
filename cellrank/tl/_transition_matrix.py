"""Transition matrix module."""

from typing import TypeVar, Iterable, Optional

from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl.kernels._velocity_kernel import BackwardMode, VelocityMode
from cellrank.tl.kernels._velocity_schemes import Scheme

AnnData = TypeVar("AnnData")


@inject_docs(m=VelocityMode, b=BackwardMode, s=Scheme)  # don't swap the order
@d.dedent
def transition_matrix(
    adata: AnnData,
    backward: bool = False,
    vkey: str = "velocity",
    xkey: str = "Ms",
    conn_key: str = "connectivities",
    gene_subset: Optional[Iterable] = None,
    mode: str = VelocityMode.DETERMINISTIC.s,
    backward_mode: str = BackwardMode.TRANSPOSE.s,
    scheme: str = Scheme.CORRELATION.s,
    softmax_scale: Optional[float] = None,
    weight_connectivities: Optional[float] = 0.2,
    density_normalize: bool = True,
    key: Optional[str] = None,
    **kwargs,
) -> KernelExpression:
    """
    Compute a transition matrix based on a combination of RNA Velocity and transcriptomic or spatial similarity.

    To learn more about the way in which the transition matrices are computed, see
    :class:`cellrank.tl.kernels.VelocityKernel` for the velocity-based transition matrix and
    :class:`cellrank.tl.kernels.ConnectivityKernel` for the similarity-based transition matrix.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    vkey
        Key from ``adata.layers`` to access the velocities.
    xkey
        Key in ``adata.layers`` where expected gene expression counts are stored.
    conn_key
        Key in :attr:`anndata.AnnData.obsp` to obtain the connectivity matrix, describing cell-cell similarity. Only
        used when ``weight_connectivities > 0``.
    gene_subset
        List of genes to be used to compute transition probabilities.
        By default, genes from ``adata.var['velocity_genes']`` are used.
    %(velocity_mode)s
    %(velocity_backward_mode_high_lvl)s
    %(velocity_scheme)s
    %(softmax_scale)s
    weight_connectivities
        Weight given to similarities as opposed to velocities. Must be in `[0, 1]`.
    density_normalize
        Whether to use density correction when computing the transition probabilities based on similarities.
        Density correction is done as by [Haghverdi16]_.
    %(write_to_adata.parameters)s
    kwargs
        Keyword arguments for :meth:`cellrank.tl.kernels.VelocityKernel.compute_transition_matrix`.

    Returns
    -------
    :class:`cellrank.tl.KernelExpression`
        A kernel expression object containing the computed transition matrix.

        %(write_to_adata)s
    """

    def compute_velocity_kernel() -> VelocityKernel:
        return VelocityKernel(
            adata, backward=backward, vkey=vkey, xkey=xkey, gene_subset=gene_subset
        ).compute_transition_matrix(
            softmax_scale=softmax_scale,
            mode=mode,
            backward_mode=backward_mode,
            scheme=scheme,
            **kwargs,
        )

    if weight_connectivities is not None:
        if 0 < weight_connectivities < 1:
            vk = compute_velocity_kernel()
            logg.info(
                f"Using a connectivity kernel with weight `{weight_connectivities}`"
            )
            ck = ConnectivityKernel(
                adata, backward=backward, conn_key=conn_key
            ).compute_transition_matrix(density_normalize=density_normalize)
            final = (
                (1 - weight_connectivities) * vk + weight_connectivities * ck
            ).compute_transition_matrix()
        elif weight_connectivities == 0:
            final = compute_velocity_kernel()
        elif weight_connectivities == 1:
            final = ConnectivityKernel(
                adata,
                backward=backward,
                conn_key=conn_key,
            ).compute_transition_matrix(density_normalize=density_normalize)
        else:
            raise ValueError(
                f"Parameter `weight_connectivities` must be in range `[0, 1]`, found `{weight_connectivities}`."
            )
    else:
        final = vk

    final.write_to_adata(key=key)

    return final
