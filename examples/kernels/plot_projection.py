"""
Plot projection
---------------

This example shows how to project a transition matrix onto a low-dimensional embedding.
"""

import scvelo as scv
import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we create the kernel.
vk = cr.tl.kernels.VelocityKernel(adata).compute_transition_matrix(
    softmax_scale=4, show_progress_bar=False
)
ck = cr.tl.kernels.ConnectivityKernel(adata).compute_transition_matrix()
k = 0.8 * vk + 0.2 * ck
k

# %%
# Next, we compute the projection of the transition matrix onto an embedding, in this case UMAP.
# The projected transition matrix can be found in ``adata.obsm['T_fwd_umap']``.
k.plot_projection(basis="umap")

# %%
# The projected transition matrix can now be visualized using :mod:`scvelo`.
# Below we use :func:`scvelo.pl.velocity_embedding_stream` to visualize it as streamlines.
scv.pl.velocity_embedding_stream(adata, vkey="T_fwd", basis="umap")

# %%
# Another way of visualizing the transition matrix is to simulate random walks on Markov chain and
# plot them in an embedding, see :ref:`sphx_glr_auto_examples_kernels_plot_random_walks.py` for more information.
