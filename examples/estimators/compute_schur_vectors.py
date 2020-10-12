# -*- coding: utf-8 -*-
"""
Compute Schur vectors
---------------------

This example shows how to compute and plot the `Schur vectors <https://en.wikipedia.org/wiki/Schur_decomposition>`_.

In CellRank, we work with the real Schur decomposition because the transition matrices we are concerned with are non
symmetrical (think of the velocity flow). In general, for these matrices the eigendecomposition does not exist.
You may still compute the top eigenvectors, however, these will in general be complex valued and hard to interpret.
The Schur decomposition is the natural generalization of the eigendecomposition for non-symmetric matrices and it’s
real variant ensures that the top Schur vectors are only real and may be interpreted in the context of metastabilities.
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
# To compute the Schur vectors, simply run the code below. The parameter ``n_components`` controls how many vectors
# to compute. We can also specify prior initial distribution over all cells as ``initial_distribution``, by default
# it is uniform.
g.compute_schur(n_components=4)

# %%
# Lastly, we plot the vectors in the UMAP embedding (default).
g.plot_schur()

# %%
# Note that above, only 3 vectors are shown, because the 1st Schur vector is a constant vector of ones.
# This can be verified visually by plotting the first Schur vector.
g.plot_schur(use=[0])

# %%
# Method :meth:`cellrank.tl.estimators.GPCCA.compute_schur` also computes the Schur matrix,
# see :ref:`sphx_glr_auto_examples_estimators_compute_schur_matrix.py`.
