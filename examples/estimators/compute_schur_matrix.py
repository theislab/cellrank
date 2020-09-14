# -*- coding: utf-8 -*-
"""
Compute Schur matrix
--------------------

This examples show how to compute and plot the Schur matrix.
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
# To compute the Schur vectors, run the code below. Parameter ``n_components`` controls how many vectors to compute.
g.compute_schur(n_components=6)

# %%
# Finally, we are ready to plot the Schur matrix. The real Schur matrix is quasi-upper triangular,
# which means there may be 1x1 or 2x2 blocks on the diagonal.
g.plot_schur_matrix()

# %%
# For more information about the Schur vectors, see :ref:`sphx_glr_auto_examples_estimators_compute_schur_vectors.py`.
