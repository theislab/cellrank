# -*- coding: utf-8 -*-
"""
Plot heatmap
------------

This example show how to plot smoothed gene expression using a heatmap.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the absorption probabilities and select a model that will be used for gene trend smoothing.
cr.tl.terminal_states(
    adata,
    cluster_key="clusters",
    weight_connectivities=0.2,
    n_states=3,
    softmax_scale=4,
    show_progress_bar=False,
)
cr.tl.lineages(adata)

model = cr.ul.models.GAM(adata)

# %%
# We can now plot the heatmap. By default, the model is fitted for all specified genes and all lineages. Paramterer
# ``show_absorption_probabilities`` is used to create a bar on top of the heatmap.
cr.pl.heatmap(
    adata,
    model,
    adata.var_names[:15],
    time_key="dpt_pseudotime",
    lineages="Alpha",
    show_absorption_probabilities=True,
    show_progress_bar=False,
)

# %%
# Sometimes, it might be beneficial to compare the smoothed expression across lineages. Parameter ``keep_gene_order``
# keeps the genes in the order as defined by the order in the first heatmap, which is the first listed lineage.
#
# Apart from the default gene sorting, we can use hierarchical clustering to cluster the genes by specifying
# ``cluster_genes=True``. We can also return the sorted/clustered genes by specifying ``return_genes=True``.
genes = cr.pl.heatmap(
    adata,
    model,
    adata.var_names[:15],
    time_key="dpt_pseudotime",
    keep_gene_order=True,
    lineages=["Alpha", "Beta"],
    show_absorption_probabilities=True,
    show_progress_bar=False,
)
genes

# %%
# Finally, we plot a heatmap-like plot where we group by genes, instead of lineages. In the case below,
# we also don't scale the expression to 0-1 range, which is the default behaviour as seen above.
cr.pl.heatmap(
    adata,
    model,
    adata.var_names[:3],
    mode="genes",
    time_key="dpt_pseudotime",
    scale=False,
    show_progress_bar=False,
)
