"""
Plot graph structures
---------------------

This functions show how to plot graph structures, such as the transition matrix.
"""

import numpy as np

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata


# %%
# First, we create a forward transition matrix using the high-level pipeline.
cr.tl.transition_matrix(
    adata, show_progress_bar=False, weight_connectivities=0.2, softmax_scale=4
)

# %%
# We can now plot the transition matrix. Below we don't show any arrows, which dramatically speeds up the plotting.
cr.pl.graph(
    adata,
    "T_fwd",
    edge_alpha=0.1,
    node_size=5,
    arrows=False,
    keys="clusters",
    keylocs="obs",
)

# %%
# To further illustrate the functionalities, let us only consider the `'Delta`' cluster. We can also filter the edges
# by their weights, as shown below. Only transitions with probability at least 0.1 are plotted.
ixs = np.where(adata.obs["clusters"] == "Delta")[0]
cr.pl.graph(adata, "T_fwd", ixs=ixs, arrows=True, node_size=200, filter_edges=(0.1, 1))

# %%
# Lastly, we can visualize different edge aggregations, such as minimum or maximum. Here we take at most 3 outgoing
# edges restricted to ``ixs`` for each node in descending order and color the nodes by the maximum outgoing weights.
# Aggregated values are always computed before any filtering happens, such as shown above.
#
# Here we also specify ``edge_reductions_restrict_to_ixs`` (by default, it is the same as ``ixs``) that computes the
# statistic between the cells marked with ``ixs`` and ``edge_reductions_restrict_to_ixs``.
#
# Below we compare the maximum transition from each of the `"Delta"` cells to any of the `"Beta"` cells.
cr.pl.graph(
    adata,
    "T_fwd",
    ixs=ixs,
    edge_alpha=0.5,
    node_size=200,
    keys="outgoing",
    arrows=False,
    top_n_edges=(3, False, "outgoing"),
    title="outgoing to Beta",
    edge_reductions=np.max,
    edge_reductions_restrict_to_ixs=np.where(adata.obs["clusters"] == "Beta")[0],
)
