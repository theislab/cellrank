"""
Compute Schur matrix
--------------------

This examples show how to compute and plot the `Schur matrix <https://en.wikipedia.org/wiki/Schur_decomposition>`_.

The Schur matrix is equivalent to the diagonal matrix of eigenvalues for diagonalizable matrices. For the real
Schur decomposition, it is a real, quasi-upper triangular matrix with 1x1 and 2x2 blocks on the diagonal.
For a rectangular, real matrix :math:`T`, it holds: :math:`T = Q R Q^T`, where :math:`R` is the real Schur matrix and
:math:`Q` is an orthogonal matrix containing the real Schur vectors as columns. 1x1 blocks on the diagonal of :math:`R`
correspond to real eigenvalues of :math:`T` whereas 2x2 blocks on the diagonal of :math:`R` correspond to pairs of
complex conjugate eigenvalues of :math:`T`.
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
# To compute the Schur vectors, run the code below. Parameter ``n_components`` controls how many vectors to compute.
g.compute_schur(n_components=6)

# %%
# Finally, we are ready to plot the Schur matrix. The real Schur matrix is quasi-upper triangular,
# which means there may be 1x1 or 2x2 blocks on the diagonal.
g.plot_schur_matrix()

# %%
# For more information about the Schur vectors, see :ref:`sphx_glr_auto_examples_estimators_compute_schur_vectors.py`.
