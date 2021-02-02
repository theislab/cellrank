"""
Compute eigendecomposition
--------------------------

This example shows how to compute and plot the spectrum of a transition matrix.

The spectrum of a rectangular matrix :math:`T` is given by its `eigenvalues
<https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix>`_. For a general, real valued matrix :math:`T`,
eigenvalues will either be real or appear in pairs of complex conjugates. CellRank can plot the spectrum in the complex
plane or restrict the visualization to the real part of the eigenvalues, which may be helpful to spot the
`eigengap <https://en.wikipedia.org/wiki/Eigengap>`_.
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
# In this example, we visualize the eigenvectors in an embedding, therefore we specify ``only_evals=False``,
# but for the method :meth:`cellrank.tl.estimators.GPCCA.plot_spectrum`, they are not necessary, because
# there we plot only the eigenvalues.
g.compute_eigendecomposition(k=20, only_evals=False)
g.eigendecomposition

# %%
# Let's start first by plotting the real spectrum. The `eigengap` is shown by the vertical line. Below
# we plot only the first 10 real eigenvalues, sorted by their real part.
g.plot_spectrum(10, real_only=True)

# %%
# Next, we plot the eigenvalues in the complex plane.
g.plot_spectrum(real_only=False)

# %%
# Finally, we can plot the left and right eigenvectors in an embedding. We can also specify which or how many
# eigenvectors to plot using the parameter ``use``. If not specified, vectors up to the `eigengap` are shown.
g.plot_eigendecomposition(left=False)
g.plot_eigendecomposition(left=True, use=2)

# %%
# Left and right eigenvectors are used in :class:`cellrank.tl.estimators.CFLARE` estimator to compute the initial or
# the states of the process.
