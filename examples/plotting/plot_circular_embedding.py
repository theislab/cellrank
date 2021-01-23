"""
Plot circular embedding
-----------------------

This example shows how to plot absorption probabilities using circular a posteriori projection as used in [Velten17]_.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the absorption probabilities.
cr.tl.terminal_states(
    adata,
    cluster_key="clusters",
    weight_connectivities=0.2,
    n_states=3,
    softmax_scale=4,
    show_progress_bar=False,
)
cr.tl.lineages(adata)

# %%
# We can now visualize the absorption probabilities by projecting them onto a unit circle. The tips of the simplex
# indicate the probability of 1.0 for the lineages and the midpoints of the edges of the edges mark where the
# probabilities of the lineages connected by an edge are equal.
cr.pl.circular_projection(
    adata, keys=["clusters", "to_terminal_states_dp"], legend_loc="upper right"
)
