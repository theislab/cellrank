# -*- coding: utf-8 -*-
"""Transition matrix module."""

from typing import Optional

from scanpy import logging as logg
from anndata import AnnData

from cellrank.tools.kernels._kernel import (
    VelocityKernel,
    KernelExpression,
    ConnectivityKernel,
)


def transition_matrix(
    adata: AnnData,
    vkey: str = "velocity",
    backward: bool = False,
    weight_connectivities: Optional[float] = None,
    sigma_corr: Optional[float] = None,
    scale_by_variances: bool = False,
    var_key: Optional[str] = "velocity_graph_uncertainties",
    var_min: float = 0.1,
    use_negative_cosines: bool = True,
    self_transitions: bool = False,
    perc: Optional[float] = None,
    threshold: Optional[float] = None,
    density_normalize: bool = True,
) -> KernelExpression:
    """
    Compute a transition matrix based on a combination of RNA Velocity and transcriptomic similarity.

    To learn more about the way in which the transition matrices are computed, see
    :class:`cellrank.tl.kernels.VelocityKernel` for the velocity-based transition matrix and
    :class:`cellrank.tl.kernels.ConnectivityKernel` for the transcriptomic-similarity-based transition matrix.

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
    use_negative_cosines
        Whether to use correlations with cells that have an angle > 90 degree with :math:`v_i`.
    sigma_corr
        Scaling parameter for the softmax. Larger values will lead to a more concentrated distribution (more peaked).
        Default is to use `1 / median_velocity_correlation`.
    scale_by_variances
        Use velocity variances to scale the softmax.
    var_key
        Key from `adata.uns` to access velocity variances.
    var_min
        Variances are clipped to this value at the lower end.
    self_transitions
        Assigns elements to the diagonal of the velocity-graph based on a confidence measure
    perc
        Quantile of the distribution of exponentiated velocity correlations. This is used as a threshold to set
        smaller values to zero.
    threshold
        Set a threshold to remove exponentiated velocity correlations smaller than :paramref:`threshold`.
    density_normalize
        Whether to use density correction when computing the transition probabilities.
        Density correction is done as by [Haghverdi16]_.

    Returns
    -------
    :class:`cellrank.tl.KernelExpression`
        A kernel expression object.
    """

    # initialise the velocity kernel and compute transition matrix
    vk = VelocityKernel(
        adata,
        backward=backward,
        vkey=vkey,
        use_negative_cosines=use_negative_cosines,
        var_key=var_key,
    )
    vk.compute_transition_matrix(
        sigma_corr=sigma_corr,
        scale_by_variances=scale_by_variances,
        var_min=var_min,
        self_transitions=self_transitions,
        perc=perc,
        threshold=threshold,
        density_normalize=density_normalize,
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
