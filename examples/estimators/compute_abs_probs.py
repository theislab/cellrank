# -*- coding: utf-8 -*-
"""
Compute absorption probabilities
--------------------------------

This example shows how to compute and plot absorption probabilities and time to absorption.
"""

import cellrank as cr
import scvelo as scv

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we prepare the kernel using high-level pipeline and the :class:`cellrank.tl.estimators.GPCCA` estimator.
k = cr.tl.transition_matrix(
    adata, weight_connectivities=0.2, softmax_scale=4, show_progress_bar=False
)
g = cr.tl.estimators.GPCCA(k)

# %%
# We need to compute or set the final states. In detail guide for both of our estimators can be found here
# :ref:`sphx_glr_auto_examples_estimators_compute_final_states_gpcca.py` or here
# :ref:`sphx_glr_auto_examples_estimators_compute_final_states_cflare.py`
g.compute_schur(n_components=4)
g.compute_metastable_states(cluster_key="clusters")
g.set_final_states_from_metastable_states(["Alpha", "Beta", "Epsilon"])

# %%
# :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities` easily scales to 100K+ cells,
# thanks to the linear solvers from :mod:`PETSc`. If needed, we can compute the absorption probabilities only to a
# few selected final states. Here, we compute the probabilities of being absorbed in all the final states.
g.compute_absorption_probabilities()

# %%
# The absorption probabilities can be inspected as seen below. Curious reader is encouraged to take a look at
# some niche tricks for :class:`cellrank.tl.Lineage` in :ref:`sphx_glr_auto_examples_other_compute_lineage_tricks.py`.
g.absorption_probabilities

# %%
# We can now plot the absorption probabilities. We can use parameters like ``same_plot`` or ``discrete`` to control
# whether to plot each lineage in a separate plot or show only the top ``n_cells`` - this is set when computing/setting
# the final states. By default, 30 cells are selected from each final state.
g.plot_absorption_probabilities()

# %%
# :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities` can also be used to compute the mean
# and the variance of time to absorption to all or just to a subset of final states.
#
# This can be specified by supplying ``time_to_absorption`` parameter. Below we compute only the mean time to
# absorption to all final states. To compute the mean and the variance only for the `'Alpha'` absorbing state,
# one specify the following ``time_to_absorption={"Alpha": "var"}``.
g.compute_absorption_probabilities(time_to_absorption="all")
g.lineage_absorption_times

# %%
# Lastly, we plot the above computed times.
adata.obs["mean_time_to_absorption"] = g.lineage_absorption_times[
    "Alpha, Beta, Epsilon mean"
]
scv.pl.scatter(adata, color="mean_time_to_absorption")
