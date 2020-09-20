# -*- coding: utf-8 -*-
"""
Plot graph structures
---------------------

This functions show how to plot graph structures, such as the transition matrix.
"""

import cellrank as cr
import numpy as np

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata


# %%
# First, we create the forward transition matrix using high-level pipeline.
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
    show_arrows=False,
    keys=("clusters",),
    keylocs="obs",
)

# %%
# To further illustrate the functionalities, let us only consider the `'Delta`' cluster. We can also filter the edges
# by their weights, as shown below. Only transition with probability at least 0.1 are plotted.
ixs = np.where(adata.obs["clusters"] == "Delta")[0]
cr.pl.graph(
    adata, "T_fwd", ixs=ixs, show_arrows=True, node_size=200, filter_edges=(0.1, 1)
)

# %%
# Below we plot at most the top 3 incoming edges restricted to ``ixs`` in descending order. Nodes are also colored in
# by the sum of incoming weights, in this case, the transition probabilities. Note that the aggregate values are
# calculated before any filtering or normalization takes place.
cr.pl.graph(
    adata,
    "T_fwd",
    ixs=ixs,
    node_size=200,
    show_arrows=True,
    top_n_edges=(3, False, "incoming"),
)

# %%
# Lastly, we can visualize different edge aggregations, such as minimum or maximum. Here we take at most 5 outgoing
# edges restricted to ``ixs`` for each node in descending order and color the nodes by the maximum outgoing weights.
#
# Here we also specify ``edge_reductions_restrict_to_ixs`` (by default, it is the same as ``ixs``) that computes the
# statistic between the cells marked with ``ixs`` and ``edge_reduction_indices``.
#
# Below we compare the maximum transition from each of the `"Delta"` cells to any of the `"Alpha"` cells or the `"Beta"`
# cells. We can visually inspect that the values are larger when going to the `"Beta"` cells, which is consistent with
# known biology, as well as our prediction, see .
cr.pl.graph(
    adata,
    "T_fwd",
    ixs=ixs,
    edge_alpha=0.5,
    node_size=200,
    keys="outgoing",
    show_arrows=False,
    top_n_edges=(5, False, "outgoing"),
    title="outgoing to Beta",
    edge_reductions=np.max,
    edge_reductions_restrict_to_ixs=np.where(adata.obs["clusters"] == "Beta")[0],
)

cr.pl.graph(
    adata,
    "T_fwd",
    ixs=ixs,
    edge_alpha=0.5,
    node_size=200,
    keys="outgoing",
    show_arrows=False,
    top_n_edges=(5, False, "outgoing"),
    title="outgoing to Alpha",
    edge_reductions=np.max,
    edge_reductions_restrict_to_ixs=np.where(adata.obs["clusters"] == "Alpha")[0],
)
