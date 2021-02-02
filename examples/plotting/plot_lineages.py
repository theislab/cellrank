"""
Plot lineages
-------------

This example shows how to plot the absorption probabilities of the process.

CellRank computes absorption probabilities to estimate how likely each individual cell is to transition into each of the
identified terminal or intermediate states.
`Absorption probabilities <https://en.wikipedia.org/wiki/Absorbing_Markov_chain>`_ refer to the probability of a random
walk starting in cell :math:`i` to reach terminal state :math:`j` before reaching any other terminal state.

Throughout CellRank, we use the terms lineage probabilities, fate probabilities and absorption probabilities
interchangeably to describe the same set of probabilities.
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
# We can now plot the terminal states. By default, we plot the degree of membership, which is available only to the
# :class:`cellrank.tl.estimators.GPCCA` estimator.
cr.pl.lineages(adata)

# %%
# We can also plot only a subset of these lineages or plot the most likely cells.
cr.pl.lineages(adata, ["Alpha", "Beta"], discrete=True)

# %%
# Lastly, we can also plot the absorption probabilities separately, one plot for each lineage.
# By default this also shows the differentiation potential, defined in [Setty19]_ as the entropy over
# the absorption probabilities.
cr.pl.lineages(adata, same_plot=False)

# %%
# To see how to compute and plot the lineage driver genes,
# refer to :ref:`sphx_glr_auto_examples_plotting_plot_lineage_drivers.py`.
