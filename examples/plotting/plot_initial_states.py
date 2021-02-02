"""
Plot initial states
-------------------

This example shows how to compute and plot the initial states of the cell-state transition.

CellRank can be applied to any cell-state transition, be it differentiation, regeneration, reprogramming or other.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the initial states. By default, we're using the :class:`cellrank.tl.estimators.GPCCA` estimator.
# The parameter ``cluster_key`` tries to associate the names of the initial states with cluster labels, whereas
# ``n_cells`` controls how many cells we take from each initial state as categorical observation - this is only
# available to the above mentioned estimator. We can show some plots of interest by specifying ``show_plots=True``.
cr.tl.initial_states(
    adata,
    cluster_key="clusters",
    n_cells=30,
    softmax_scale=4,
    n_states=1,
    show_progress_bar=False,
)

# %%
# We can now plot the initial states. By default, when using :class:`cellrank.tl.estimators.GPCCA`, we plot continuous
# membership vectors to visualize individual cells associations with an initial state.
#
# We can also plot membership vectors for different initial states separately if we computed more than one initial
# state using ``same_plot=False``. As we only have one initial state here, this does not make sense.
cr.pl.initial_states(adata)

# %%
# Lastly, we can discretize the assignment of cells to initial states by showing the cells most
# likely to belong to the initial state by specifying the ``discrete`` parameter.
cr.pl.initial_states(adata, discrete=True)

# %%
# To see how to compute and plot the terminal states or the lineages, see
# :ref:`sphx_glr_auto_examples_plotting_plot_terminal_states.py` or
# :ref:`sphx_glr_auto_examples_plotting_plot_lineages.py`, respectively.
