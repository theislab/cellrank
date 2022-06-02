"""
Compute coarse-grained transition matrix
----------------------------------------

This example shows how to compute and plot the coarse-grained transition matrix among the set of macrostates.

We computed the macrostates using Generalized Perron Cluster Cluster Analysis :cite:`reuter:18` :cite:`reuter:19`,
see the example :ref:`sphx_glr_auto_examples_estimators_compute_macrostates.py`. The coarse-grained transition matrix
shows transitional behavior among the set of macrostates and may be used to classify these states as initial,
intermediate or terminal.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we prepare the kernel and the :class:`cellrank.tl.estimators.GPCCA` estimator.
vk = cr.tl.kernels.VelocityKernel(adata).compute_transition_matrix(
    softmax_scale=4, show_progress_bar=False
)
ck = cr.tl.kernels.ConnectivityKernel(adata).compute_transition_matrix()
k = 0.8 * vk + 0.2 * ck
g = cr.tl.estimators.GPCCA(k)

# %%
# Next, we compute the Schur vectors.
g.compute_schur(n_components=6)

# %%
# Now we can compute the macrostates of the process, which also computes the coarse-grained transition matrix.
# Here, parameter ``cluster_key`` tries to associate the names of macrostates with the cluster labels.
# For more options, see :ref:`sphx_glr_auto_examples_estimators_compute_macrostates.py`.
g.compute_macrostates(n_states=6, cluster_key="clusters")

# %%
# We can now plot the coarse-grained transition matrix.
g.plot_coarse_T(text_kwargs={"fontsize": 10})

# %%
# The coarse-grained transition matrix can also be used when setting the terminal states, see
# :ref:`sphx_glr_auto_examples_estimators_compute_terminal_states_gpcca.py`.
