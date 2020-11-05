# -*- coding: utf-8 -*-
"""
Compute potential lineage drivers
---------------------------------

This example shows how to compute and plot expression trends for genes which may be involved in lineage decisions.

We identify these by correlating gene expression with absorption probabilities towards a specific terminal state.
"""

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
# We need to compute the absorption probabilities. In this example, we're using
# :class:`cellrank.tl.estimators.GPCCA` estimator to estimate the terminal states of the process, but
# :class:`cellrank.tl.estimators.CFLARE` can be used as well.
#
# In detail guide for :class:`cellrank.tl.estimators.GPCCA` estimator can be found in
# :ref:`sphx_glr_auto_examples_estimators_compute_terminal_states_gpcca.py`.
g.compute_schur(n_components=4)
g.compute_macrostates(cluster_key="clusters")
g.set_terminal_states_from_macrostates(["Alpha", "Beta", "Epsilon"])
g.compute_absorption_probabilities()
g.absorption_probabilities

# %%
# To compute the potential driver genes, simply call the
# :meth:`cellrank.tl.estimators.BaseEstimator.compute_lineage_drivers` method. By default, the these are computed for
# all lineages. We can restrict this computation to only a few clusters, using the ``cluster_key`` and ``clusters``
# parameters.
# We also compute the corrected p-values (qval) and the 95% confidence intervals for the correlations.
g.compute_lineage_drivers(lineages="Alpha")
g.lineage_drivers.sort_values("Alpha corr", ascending=False)

# %%
# Lastly, we plot the top 3 potential driver genes for the `"Alpha"` lineage.
g.plot_lineage_drivers("Alpha", n_genes=3)
