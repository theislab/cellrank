# -*- coding: utf-8 -*-
"""
Compute macrostates
-------------------

This example shows how to compute and plot the macrostates.

For the computation of macrostates, we adapted the popular Generalized Perron Cluster Cluster Analysis [GPCCA18]_
[Reuter19]_ method to the single cell context. We provide a scalable implementation which can decompose datasets of
100k+ cells into their dominant dynamical macrostates in just a few minutes. GPCCA relies on the real Schur
decomposition to handle non-symmetric transition matrices as they arise from RNA velocity information, see
:ref:`sphx_glr_auto_examples_estimators_compute_schur_vectors.py` and
:ref:`sphx_glr_auto_examples_estimators_compute_schur_matrix.py`.
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
# First, we need to compute the Schur vectors. By default, the first 10 vectors are computed.
g.compute_schur(n_components=4)

# %%
# We can now compute the macrostates of the Markov chain. The first important parameter is the ``cluster_key``,
# which tries to associate the names of macrostates with the cluster labels.
# The second import parameter is ``n_cells``, which selects the top cells from each state based
# on the membership degree. By default, 30 cells are selected.
#
# Lastly, the parameter ``n_states`` can also be estimated by using either the `eigengap` or the `minChi` criterion from
# [GPCCA18]_.
g.compute_macrostates(n_states=3, cluster_key="clusters")

# %%
# After computing the macrostates, we can inspect them as follows. Below we show for each cell the membership
# degree of the macrostates.
g.macrostates_memberships

# %%
# To get the categorical observations, top ``n_cells`` for each macrostate, we can inspect the attribute below.
g.macrostates

# %%
# We can now plot the membership degree, as well as the categorical assignment.
g.plot_macrostates()
g.plot_macrostates(discrete=True)

# %%
# Both of these options are shown in the same plot, which is not always desirable.
# To change this, simply run the following.
g.plot_macrostates(same_plot=False)

# %%
# Method :meth:`cellrank.tl.estimators.GPCCA.compute_macrostates` also computes the coarse-grained transition
# matrix between the macrostates, see :ref:`sphx_glr_auto_examples_estimators_compute_coarse_T.py`.
