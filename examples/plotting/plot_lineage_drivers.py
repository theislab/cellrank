# -*- coding: utf-8 -*-
"""
Plot potential lineage drivers
------------------------------

This example shows how to compute and plot expression trends for genes which may be involved in lineage decisions.

We identify these by correlating gene expression with absorption probabilities towards a specific terminal state.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we need to compute the terminal states and the absorption probabilities towards them.
cr.tl.terminal_states(
    adata,
    cluster_key="clusters",
    n_cells=30,
    n_states=3,
    softmax_scale=4,
    show_progress_bar=False,
)
cr.tl.lineages(adata)

# %%
# Once the lineages have been computed, we can compute the potential driver genes for each of them. It is also
# possible to restrict this computation to just a few clusters, defined by ``cluster_key`` and ``clusters``.
#
# By default we are computing the driver genes for all lineages.
drivers = cr.tl.lineage_drivers(adata, lineages="Alpha")
drivers

# %%
# Finally, we can plot the potential drivers. Below we plot top 3 driver genes for the `'Alpha'` lineage.
cr.pl.lineage_drivers(adata, lineage="Alpha", n_genes=3)
