# -*- coding: utf-8 -*-
"""
Plot lineages
-------------

This example shows how to plot absorption probabilities of the process.
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
# We can now plot the initial states. By default, we plot the degree of membership, which is available only to the
# :class:`cellrank.tl.estimators.GPCCA` estimator.
cr.pl.lineages(adata)

# %%
# We can also plot only a subset of these lineages or plot the most likely cells.
cr.pl.lineages(adata, ["Alpha", "Beta"], discrete=True)

# %%
# Lastly, we can also plot the absorption probabilities separately, one plot for each lineage.
# By default this also shows the differentiation potential.
cr.pl.lineages(adata, same_plot=False)

# %%
# To see how to compute and plot the lineage driver genes,
# refer to :ref:`sphx_glr_auto_examples_plotting_plot_lineage_drivers.py`.
