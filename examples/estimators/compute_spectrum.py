# -*- coding: utf-8 -*-
"""
Compute spectrum
----------------

This example shows how to compute and plot the spectrum of the transition matrix.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we the kernel using high-level pipeline and the :class:`cellrank.tl.estimators.GPCCA` estimator.
k = cr.tl.transition_matrix(
    adata, weight_connectivities=0.2, softmax_scale=4, show_progress_bar=False
)
g = cr.tl.estimators.GPCCA(k)

# %%
# In this example, we iwant visualize the eigenvectors in an embedding, therefore we specify ``only_evals=False``,
# but for the method :meth:`cellrank.tl.estimators.BaseEstimator.plot_spectrum`, they are not necessary.
g.compute_eigendecomposition(k=20, only_evals=False)
g.eigendecomposition

# %%
# Let's start first by plotting the real spectrum. The `eigengap` is shown by the vertical line. Below
# we plot only the first 10 real eigenvalues, sorted by their real part.
g.plot_spectrum(10, real_only=True)

# %%
# And similarly for the complex spectrum.
g.plot_spectrum(real_only=False)

# %%
# Finally, we can plot the left and right eigenvectors in an embedding. We can also specify which or how many
# eigenvectors to plot using the parameter ``use``. If not specified, vectors up to the `eigengap` are shown.
g.plot_eigendecomposition(left=False)
g.plot_eigendecomposition(left=True, use=2)

# %%
# Left and right eigenvectors are used in :class:`cellrank.tl.estimators.CFLARE` estimator to compute the final states
# of the process. For more information, see :ref:`sphx_glr_auto_examples_estimators_compute_final_states_cflare.py`.
