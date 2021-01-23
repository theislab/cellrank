"""
Compute absorption probabilities
--------------------------------

This example shows how to compute and plot the absorption probabilities and the time to absorption.

Absorption probabilities are used in CellRank to define how likely each individual cell is to transition into each of
the identified terminal states. In a usual workflow, you would first compute a set of terminal states for your dataset
and next ask the question how likely each cell is to develop towards each of these states. CellRank provides an
efficient implementation of computing the absorption probabilities that scales to 100k+ cells.
"""

import scvelo as scv

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we prepare the kernel using the high-level pipeline and the :class:`cellrank.tl.estimators.GPCCA` estimator.
k = cr.tl.transition_matrix(
    adata, weight_connectivities=0.2, softmax_scale=4, show_progress_bar=False
)
g = cr.tl.estimators.GPCCA(k)

# %%
# We need to compute or set the terminal states. In detail guide for both of our estimators can be found here
# :ref:`sphx_glr_auto_examples_estimators_compute_terminal_states_gpcca.py`.
g.compute_schur(n_components=4)
g.compute_macrostates(cluster_key="clusters")
g.set_terminal_states_from_macrostates(["Alpha", "Beta", "Epsilon"])

# %%
# :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities` easily scales to 100k+ cells,
# thanks to the linear solvers from :mod:`PETSc`.
#
# The computation of absorption probabilities may be restricted to a subset of the identified states via the ``keys``
# parameter. In our case, we are interested in the absorption probabilities towards each of the terminal states we
# identified in our data.
g.compute_absorption_probabilities()

# %%
# The absorption probabilities can be inspected as seen below. Curious reader is encouraged to take a look at
# some niche tricks for :class:`cellrank.tl.Lineage` in :ref:`sphx_glr_auto_examples_other_compute_lineage_tricks.py`.
g.absorption_probabilities

# %%
# We can now plot the absorption probabilities. We can use parameters like ``same_plot`` or ``discrete`` to control
# whether to plot each lineage in a separate plot or show only the top ``n_cells`` - this is set when computing/setting
# the terminal states. By default, 30 cells are selected from each terminal state.
g.plot_absorption_probabilities()

# %%
# :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities` can also be used to compute the mean
# and the variance of time to absorption to all or just to a subset of terminal states.
#
# This can be specified by supplying ``time_to_absorption`` parameter. Below we compute only the mean time to
# absorption to all terminal states. To compute the mean and the variance only for the `"Alpha"` absorbing state,
# one specify the following ``time_to_absorption={"Alpha": "var"}``.
g.compute_absorption_probabilities(time_to_absorption="all")
g.lineage_absorption_times

# %%
# Lastly, we plot the above computed time.
adata.obs["mean_time_to_absorption"] = g.lineage_absorption_times[
    "Alpha, Beta, Epsilon mean"
]
scv.pl.scatter(adata, color="mean_time_to_absorption")
