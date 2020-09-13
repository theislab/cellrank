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
# First, let us prepare the kernel using high-level pipeline and the :class:`cellrank.tl.estimators.GPCCA` estimator.
k = cr.tl.transition_matrix(adata, show_progress_bar=False)
g = cr.tl.estimators.GPCCA(k)

# %%
# To compute the Schur vectors, simply run the code below. Parameter ``n_components``
# controls how many vectors to compute.
g.compute_schur(n_components=20)

# %%
# Finally, we are ready to plot the Schur matrix. The real Schur matrix is quasi-upper triangular,
# which means there may be 1x1 or 2x2 blocks on the diagonal.
g.plot_schur_matrix()
