# -*- coding: utf-8 -*-
"""
Plot initial states
-------------------

This example shows how to compute and plot the initial states of the process.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the initial states. By default, we're using the :class:`cellrank.tl.estimators.GPCCA` estimator.
# Parameter ``cluster_key`` tries to associate the names of the initial states with cluster labels, whereas
# ``n_cells`` controls how many cells we take from each final state as categorical observation - this is only available
# for the above mentioned estimator. Plots showing the terminal states and more can be specified by ``show_plots=True``.
cr.tl.initial_states(
    adata,
    cluster_key="clusters",
    n_cells=30,
    softmax_scale=4,
    n_states=1,
    show_progress_bar=False,
)

# %%
# We can now plot the initial states. By default, we plot the membership degree, which is only available to
# :class:`cellrank.tl.estimators.GPCCA` estimator.
cr.pl.initial_states(adata)

# %%
# Lastly, we can also plot the discrete values, as seen below.
cr.pl.initial_states(adata, discrete=True)

# %%
# To see how to compute and plot the terminal states or the lineages, see
# :ref:`sphx_glr_auto_examples_plotting_plot_terminal_states.py` and
# see :ref:`sphx_glr_auto_examples_plotting_plot_lineages.py`.
