# -*- coding: utf-8 -*-
"""
Plot Schur vectors
------------------

This examples show how to compute and plot the Schur vectors.
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
g.compute_schur(n_components=6)

# %%
# Lastly, we plot the vectors in the UMAP embedding (default).
g.plot_schur()

# %%
# Note that above, only 5 vectors are shown, because the 1st Schur vector is a unit vector. This can be verified
# visually by plotting the first Schur vector.
g.plot_schur(use=[0])
