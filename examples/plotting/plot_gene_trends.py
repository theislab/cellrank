# -*- coding: utf-8 -*-
"""
Plot gene trends
----------------

This example shows how to plot smoothed gene expression across lineages.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the absorption probabilities and the model that will be used for gene trend smoothing.
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
# To plot the trends for some genes, run the code below. Parameter ``data_key`` specifies layer in ``adata.layers``
# from which we take our gene expression.
cr.pl.gene_trends(
    adata,
    model,
    adata.var_names[:2],
    data_key="Ms",
    show_progres_bar=False,
    time_key="dpt_pseudotime",
)

# %%
# We can plot all the trends in the same plot and hide the cells, as seen below.
cr.pl.gene_trends(
    adata,
    model,
    adata.var_names[:2],
    data_key="Ms",
    same_plot=True,
    hide_cells=True,
    time_key="dpt_pseudotime",
    show_progres_bar=False,
)

# %%
# We can specify specific pseudotime ranges, separately for each lineage - the model is always fitted
# using cells with a given time range. By default, the start is always 0 and the end is automatically determined
# from the absorption probabilities.
cr.pl.gene_trends(
    adata,
    model,
    adata.var_names[:2],
    data_key="Ms",
    lineages=["Alpha", "Beta"],
    time_range=[(0.2, 1), (0, 0.8)],
    time_key="dpt_pseudotime",
    show_progres_bar=False,
)

# %%
# There are many more options as how to customize the plot or how to pass additional arguments to
# :meth:`cellrank.ul.models.BaseModel` and we encourage the reader to take a closer look at the documentation of
# :func:`cellrank.pl.gene_trends`.
