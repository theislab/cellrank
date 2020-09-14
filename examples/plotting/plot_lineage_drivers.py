# -*- coding: utf-8 -*-
"""
Plot lineages drivers
---------------------

This example shows how to compute and plot lineage driver genes.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we need to compute the initial or terminal states and the absorption probabilities towards them.
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
# Once the lineages have been computed, we can compute the driver genes for each of them. It is also possible to
# restrict this computation to just a few clusters, defined by ``cluster_key`` and ``clusters``. By default
# we are computing the driver genes for all lineages.
cr.tl.lineage_drivers(adata)

# %%
# Finally, we can plot the driver genes. Below we plot top 3 driver genes for the `'Alpha'` lineage.
cr.pl.lineage_drivers(adata, lineage="Alpha", n_genes=3)
