"""
Plot random walks
-----------------

This example shows how to plot random walks on a Markov chain in an embedding.
"""

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
# Next, we simulate and plot random walks. Below we simulate 100 random walks of length 100 where starting points
# are uniformly sampled from the `Ngn3 low EP` cluster.
k.plot_random_walks(
    100,
    start_ixs={"clusters": "Ngn3 low EP"},
    max_iter=100,
    show_progress_bar=False,
    ixs_legend_loc="best",
    seed=42,
)


# %%
# Apart from ``start_ixs``, we can also specify ``stop_ixs`` and ``successive_hits`` to control when random walks
# stop. Random walk is stopped if the maximum number of iterations is reached or when states in ``stop_ixs``
# is visited successively ``successive_hits`` times.
k.plot_random_walks(
    100,
    start_ixs={"clusters": "Ngn3 low EP"},
    stop_ixs={"clusters": ["Alpha", "Beta"]},
    max_iter=100,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",
    seed=42,
)

# %%
# Another way of visualizing the transition matrix is to project it onto an embedding,
# see :ref:`sphx_glr_auto_examples_kernels_plot_projection.py` for more information.
