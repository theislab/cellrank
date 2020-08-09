# -*- coding: utf-8 -*-
"""Transition matrix module."""

from typing import TypeVar, Iterable, Optional

from cellrank import logging as logg
from cellrank.utils._docs import d, inject_docs
from cellrank.tools.kernels import VelocityKernel, ConnectivityKernel
from cellrank.tools.kernels._base_kernel import KernelExpression
from cellrank.tools.kernels._velocity_kernel import BackwardMode, VelocityMode

AnnData = TypeVar("AnnData")


@inject_docs(m=VelocityMode, b=BackwardMode)  # don't swap the order
@d.dedent
def transition_matrix(
    adata: AnnData,
    backward: bool = False,
    vkey: str = "velocity",
    xkey: str = "Ms",
    gene_subset: Optional[Iterable] = None,
    mode: str = VelocityMode.DETERMINISTIC.s,
    backward_mode: str = BackwardMode.TRANSPOSE.s,
    seed: Optional[int] = None,
    softmax_scale: int = 4.0,
    weight_connectivities: Optional[float] = None,
    density_normalize: bool = True,
) -> KernelExpression:
    """
    Compute a transition matrix based on a combination of RNA Velocity and transcriptomic similarity.

    To learn more about the way in which the transition matrices are computed, see
    :class:`cellrank.tl.kernels.VelocityKernel` for the velocity-based transition matrix and
    :class:`cellrank.tl.kernels.ConnectivityKernel` for the transcriptomic-similarity-based transition matrix.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    vkey
        Key from :paramref:`adata` `.layers` to access the velocities.
    xkey
        Key in :paramref:`adata` `.layers` where expected gene expression counts are stored.
    gene_subset
        List of genes to be used to compute transition probabilities. By default, the `velocity_genes` of
        :paramref:`adata` `. var` are used.
    %(velocity_mode)s
    %(velocity_backward_mode)s
    seed
        Set the seed for random state, only used when sampling.
    softmax_scale
        Scaling parameter for the softmax.
    weight_connectivities
        Weight given to transcriptomic similarities as opposed to velocities. Must be in `[0, 1]`.
    density_normalize
        Whether to use density correction when computing the transition probabilities based on connectivities.
        Density correction is done as by [Haghverdi16]_.

    Returns
    -------
    :class:`cellrank.tl.KernelExpression`
        A kernel expression object.
    """

    # initialise the velocity kernel and compute transition matrix
    vk = VelocityKernel(
        adata, backward=backward, vkey=vkey, xkey=xkey, gene_subset=gene_subset
    )
    vk.compute_transition_matrix(
        softmax_scale=softmax_scale, mode=mode, backward_mode=backward_mode, seed=seed
    )

    if weight_connectivities is not None:
        if 0 < weight_connectivities < 1:
            logg.info(
                f"Using a connectivity kernel with weight `{weight_connectivities}`"
            )
            ck = ConnectivityKernel(adata, backward=backward).compute_transition_matrix(
                density_normalize=density_normalize
            )
            final = (1 - weight_connectivities) * vk + weight_connectivities * ck
        elif weight_connectivities == 0:
            final = vk
        elif weight_connectivities == 1:
            final = ConnectivityKernel(
                adata, backward=backward
            ).compute_transition_matrix(density_normalize=density_normalize)
        else:
            raise ValueError(
                f"Parameter `weight_connectivities` must be in range `[0, 1]`, found `{weight_connectivities}`."
            )
    else:
        final = vk
    final.write_to_adata()

    return final
