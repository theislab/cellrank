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
# First, we create the kernel using the high-level function :func:`cellrank.tl.transition_matrix`, but any
# instance of :class:`cellrank.tl.kernels.Kernel` would do.
k = cr.tl.transition_matrix(
    adata, weight_connectivities=0.2, softmax_scale=4, show_progress_bar=False
)
k

# %%
# Next, we compute the projection of the transition matrix onto an embedding, in this case UMAP.
# The projected transition matrix can be found in ``adata.obsm['T_fwd_umap']``.
k.compute_projection(basis="umap")

# %%
# The projected transition matrix can now be visualized using :mod:`scvelo`.
# Below we use :func:`scvelo.pl.velocity_embedding_stream` to visualize it as streamlines.
scv.pl.velocity_embedding_stream(adata, vkey="T_fwd", basis="umap")

# %%
# Another way of visualizing the transition matrix is to simulate random walks on Markov chain and
# plot them in an embedding, see :ref:`sphx_glr_auto_examples_kernels_plot_random_walks.py` for more information.
