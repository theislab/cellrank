# -*- coding: utf-8 -*-
"""
Compute driver genes
--------------------

This example shows how to compute and plot lineage driver genes.
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
# We need to compute the absorption probabilities. In this example, we're using
# :class:`cellrank.tl.estimators.GPCCA` estimator to estimate the final states of the process, but
# :class:`cellrank.tl.estimators.CFLARE` can be used as well.
#
# In detail guide for both of our estimators can be found in
# :ref:`sphx_glr_auto_examples_estimators_compute_final_states_gpcca.py` or in
# :ref:`sphx_glr_auto_examples_estimators_compute_final_states_cflare.py`
g.compute_schur(n_components=4)
g.compute_metastable_states(cluster_key="clusters")
g.set_final_states_from_metastable_states(["Alpha", "Beta", "Epsilon"])
g.compute_absorption_probabilities()
g.absorption_probabilities

# %%
# To compute the driver genes, simply call the :meth:`cellrank.tl.estimators.BaseEstimator.compute_lineage_drivers`
# By default, lineage drivers are computed for all lineages. We can restrict this computation to only a few clusters,
# using ``cluster_key`` and ``clusters``.
g.compute_lineage_drivers()
g.lineage_drivers

# %%
# Lastly, we plot the top 3 driver genes for the `'Alpha'` lineage.
g.plot_lineage_drivers("Alpha", n_genes=3)
