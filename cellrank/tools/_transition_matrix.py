# -*- coding: utf-8 -*-
from typing import Optional
from anndata import AnnData
from scanpy import logging as logg

from cellrank.tools.kernels._kernel import (
    KernelExpression,
    VelocityKernel,
    ConnectivityKernel,
)


def transition_matrix(
    adata: AnnData,
    vkey: str = "velocity",
    backward: bool = False,
    weight_connectivities: Optional[float] = None,
    density_normalize: bool = True,
) -> KernelExpression:
    """
    High level function to compute a transition matrix based on a combination
    of RNA Velocity and transcriptomic similarity.

    To learn more about the way in which the transition matrices are computed, see
    :class:`cellrank.tl.kernels.VelocityKernel` for the velocity-based transition matrix and
    :class:`cellrank.tl.kernels.ConnectivityKernel` for the transcritomic-similarity-based transition matrix.

    Params
    ------
    adata: :class:`anndata.AnnData`
        Annotated data object.
    vkey
        Key from :paramref:`adata` `.layers` to access the velocities.
    backward
        Direction of the process.
    weight_connectivities
        Weight given to transcriptomic similarities as opposed to velocities. Must be in `[0, 1]`.
    density_normalize
        Whether to use density correction when computing the transition probabilities.
        Density correction is done [Haghverdi16]_.

    Returns
    -------
    :class:`cellrank.tl.KernelExpression`
        A kernel object.
    """

    # initialise the kernel objects
    vk = VelocityKernel(adata, backward=backward, vkey=vkey).compute_transition_matrix(
        density_normalize=density_normalize
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
                f"The parameter `weight_connectivities` must be in range `[0, 1]`, found `{weight_connectivities}`."
            )
    else:
        final = vk
    final.write_to_adata()

    return final
