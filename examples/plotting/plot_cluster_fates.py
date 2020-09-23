# -*- coding: utf-8 -*-
"""
Plot aggregated cellular fates
------------------------------

This example shows how to aggregate fate probabilities from the single cell level to a cluster level
and how to visualize these in various ways.
"""

import cellrank as cr
import scanpy as sc

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the terminal states and the lineages, as well as :mod:`scanpy`'s PAGA [Wolf19]_.
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
# Violin plots are helpful to visualize the distribution of fate probabilities per cluster.
# It is also possible to restrict this plot only to a subset of clusters using the ``clusters`` parameter.
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
cr.pl.cluster_fates(adata, mode="paga", legend_loc="on data", basis="umap")
