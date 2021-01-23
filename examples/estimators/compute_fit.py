"""
Fit estimator
-------------

This example shows how to fit an estimator in order to compute the lineages.

CellRank is composed of :mod:`cellrank.tl.kernels` and :mod:`cellrank.tl.estimators`. Kernels compute transition
matrices on the basis of directed data, given by i.e. RNA velocity [Manno18]_ [Bergen20]_ while estimators perform
inference making use of kernels.

Here, we show how to fit an estimator, given a kernel. An estimator may be used e.g. to identify initial and terminal
states or to compute lineage probabilities towards the terminal states, for each individual cell.
The :meth:`cellrank.tl.estimators.BaseEstimator.fit` method is a combination of individual methods to compute terminal
states and absorption probabilities.

To have maximum control over the inference performed, use these individual methods directly.
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
# to compute the terminal states and optionally the absorption probabilities by following similar steps as defined in
# :ref:`sphx_glr_auto_examples_estimators_compute_terminal_states_gpcca.py`.
g.fit(n_lineages=3, cluster_key="clusters", compute_absorption_probabilities=True)

# %%
# In order to verify that the absorption probabilities have been computed, we plot them below.
g.plot_absorption_probabilities()
