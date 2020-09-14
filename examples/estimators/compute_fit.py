# -*- coding: utf-8 -*-
"""
Fit estimator
-------------

This example shows how to fit an estimator in order to compute the lineages.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we prepare the kernel using high-level pipeline and the :class:`cellrank.tl.estimators.GPCCA` estimator.
k = cr.tl.transition_matrix(
    adata, weight_connectivities=0.2, softmax_scale=4, show_progress_bar=False
)
g = cr.tl.estimators.GPCCA(k)

# %%
# Afterwards, we simply call :meth:`cellrank.tl.estimators.BaseEstimator.fit`. It offers a quick and easy way
# to compute the final states and optionally also the absorption probabilities by following similar steps as defined in
# :ref:`sphx_glr_auto_examples_estimators_compute_final_states_gpcca.py` or in
# :ref:`sphx_glr_auto_examples_estimators_compute_final_states_cflare.py`, depending on the estimator.
g.fit(n_lineages=3, cluster_key="clusters", compute_absorption_probabilities=True)

# %%
# In order to verify that the absorption probabilities have been computed, we plot them below.
g.plot_absorption_probabilities()
