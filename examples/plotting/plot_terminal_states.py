"""
Plot terminal states
--------------------

This example shows how to compute and plot the terminal states of the cell-state transition.

CellRank can be applied to any cell-state transition, be it differentiation, regeneration, reprogramming or other.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the terminal states. By default, we're using the :class:`cellrank.tl.estimators.GPCCA` estimator.
# The parameter ``cluster_key`` tries to associate the names of the terminal states with cluster labels, whereas
# ``n_cells`` controls how many cells we take from each terminal state as categorical observation - this is only
# available to the above mentioned estimator. We can show some plots of interest by specifying ``show_plots=True``.
cr.tl.terminal_states(
    adata,
    cluster_key="clusters",
    n_cells=30,
    softmax_scale=4,
    n_states=3,
    show_progress_bar=False,
)

# %%
# We can now plot the terminal states. By default, when using :class:`cellrank.tl.estimators.GPCCA`, we plot continuous
# membership vectors to visualize individual cells associations with an terminal state.
cr.pl.terminal_states(adata)

# %%
# We can also plot membership vectors for different terminal states separately  state using ``same_plot=False``.
cr.pl.terminal_states(adata, same_plot=False)

# %%
# Lastly, we can discretize the assignment of cells to terminal states by showing the cells most
# likely to belong to the terminal state by specifying the ``discrete`` parameter.
cr.pl.terminal_states(adata, discrete=True)

# %%
# To see how to compute and plot the terminal states or the lineages, see
# :ref:`sphx_glr_auto_examples_plotting_plot_initial_states.py` or
# :ref:`sphx_glr_auto_examples_plotting_plot_lineages.py`, respectively.
