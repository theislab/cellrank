# -*- coding: utf-8 -*-
"""
Plot cluster fates
------------------

This example shows how to plot average absorption probabilities in a cluster specific way.
"""

import cellrank as cr
import scanpy as sc

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the terminal states and the lineages, as well as scanpy's PAGA.
cr.tl.terminal_states(
    adata,
    cluster_key="clusters",
    weight_connectivities=0.2,
    n_states=3,
    softmax_scale=4,
    show_progress_bar=False,
)
cr.tl.lineages(adata)

sc.tl.paga(adata, "clusters")

# %%
# We can visualize the aggregate absorption probabilities as a bar plot.
cr.pl.cluster_fates(adata, mode="bar")

# %%
# Similarly, the aggregate information can be visualized using a heatmap or a clustermap.
cr.pl.cluster_fates(adata, mode="heatmap")
cr.pl.cluster_fates(adata, mode="clustermap")

# %%
# Violin plot can be useful to visualize the distribution of absorption probabilities per cluster. It is also
# possible to restric this plot only to a subset of clusters.
cr.pl.cluster_fates(adata, mode="violin", cluster_key="clusters")

# %%
# We can also plot the PAGA graph, using an UMAP embedding while visualizing the aggregate absorption probabilities
# as pie charts.
#
# However, using :func:`scvelo.tl.paga` and :mod:`cellrank`'s initial and terminal state probabilities, we can also
# visualize a directed version of PAGA, see :ref:`sphx_glr_auto_examples_plotting_plot_directed_paga.py`.
cr.pl.cluster_fates(adata, mode="paga_pie", basis="umap", cluster_key="clusters")

# %%
# Lastly, we can visualize the absorption probabilities in PAGA graph by coloring the nodes.
cr.pl.cluster_fates(adata, mode="paga")
