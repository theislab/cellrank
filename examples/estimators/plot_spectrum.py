# -*- coding: utf-8 -*-
"""
Plot metastable states
----------------------

This examples show how to compute and plot metastable states obtained by :class:`cellrank.tl.estimators.GPCCA`.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, let us prepare the kernel using high-level pipeline and the :class:`cellrank.tl.estimators.GPCCA` estimator.
k = cr.tl.transition_matrix(adata, show_progress_bar=False)
g = cr.tl.estimators.GPCCA(k)

# %%
# We can decide whether we want to also compute eigenvectors (will be a bit slower). In this example, we also
# visualize the eigenvectors in an embedding, so we compute them as well, but for the method
# :meth:`cellrank.tl.estimators.BaseEstimator.plot_spectrum` they are not necessary.
g.compute_eigendecomposition(k=20, only_evals=False)
g.eigendecomposition

# %%
# Let's start first by plotting the real spectrum. The eigengap is shown by the horizontal line. Below
# we plot only the first 10 real eigenvalues.
g.plot_spectrum(10, real_only=True)

# %%
# And similarly for the complex spectrum.
g.plot_spectrum(real_only=False)

# %%
# Finally, we can plot the left and right eigenvectors in an embedding. We can also specify which or how many
# vectors to plot using the parameter ``use``. If not specified, vectors upto the eigengap are plotted.
g.plot_eigendecomposition(left=False)
g.plot_eigendecomposition(left=True, use=2)
