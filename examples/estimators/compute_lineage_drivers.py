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
# First, we prepare the kernel and the :class:`cellrank.estimators.GPCCA` estimator.
vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix(
    softmax_scale=4, show_progress_bar=False
)
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
k = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(k)

# %%
# We need to compute the absorption probabilities. In this example, we're using
# :class:`cellrank.estimators.GPCCA` estimator to estimate the terminal states of the process, but
# :class:`cellrank.estimators.CFLARE` can be used as well.
#
# In detail guide for :class:`cellrank.estimators.GPCCA` estimator can be found in
# :ref:`sphx_glr_auto_examples_estimators_compute_terminal_states_gpcca.py`.
g.compute_schur(n_components=4)
g.compute_macrostates(cluster_key="clusters")
g.set_terminal_states_from_macrostates(["Alpha", "Beta", "Epsilon"])
g.compute_absorption_probabilities()
g.absorption_probabilities

# %%
# To compute the potential driver genes, simply call the
# :meth:`cellrank.estimators.BaseEstimator.compute_lineage_drivers` method. By default, the these are computed for
# all lineages. We can restrict this computation to only a few clusters, using the ``cluster_key`` and ``clusters``
# parameters.
# We also compute the corrected p-values (qval) and the 95% confidence intervals for the correlations.
g.compute_lineage_drivers(lineages="Alpha")
g.lineage_drivers.sort_values("Alpha_corr", ascending=False)

# %%
# Lastly, we plot the top 3 potential driver genes for the `"Alpha"` lineage.
g.plot_lineage_drivers("Alpha", n_genes=3)
