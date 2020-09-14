# -*- coding: utf-8 -*-
"""
Compute coarse-grained transition matrix
----------------------------------------

This example shows how to compute and plot coarse-grained transition matrix of the metastable states.
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
# Next, we compute the Schur vectors. See :ref:`sphx_glr_auto_examples_estimators_compute_schur_vectors.py` for
# more information.
g.compute_schur(n_components=6)

# %%
# Now we can compute the metastable states of the process, which also computes the coarse-grained transition matrix.
# Here, parameter ``cluster_key`` tries to associate the names of metastable states with the cluster labels.
# For more options, see :ref:`sphx_glr_auto_examples_estimators_compute_metastable_states.py`.
g.compute_metastable_states(n_states=6, cluster_key="clusters")

# %%
# We can now plot the coarse-grained transition matrix. Apart from the stationary distribution, which is shown by
# default (if it exists), we can also plot the initial distribution of the metastable states.
g.plot_coarse_T(show_initial_dist=True, text_kwargs={"fontsize": 10})

# %%
# The coarse-grained transition matrix can also be used when setting the final states, see
# :ref:`sphx_glr_auto_examples_estimators_compute_final_states_gpcca.py`.
