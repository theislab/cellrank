"""
Compute terminal states using GPCCA
-----------------------------------

This example shows how to compute and plot the terminal states using the :class:`cellrank.estimators.GPCCA`.

This estimator makes use of Generalized Perron Cluster Cluster Analysis :cite:`reuter:18` :cite:`reuter:19`
as implemented in `pyGPCCA <https://pygpcca.readthedocs.io/en/latest/>`_.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we prepare the kernel and the :class:`cellrank.estimators.GPCCA` estimator.
vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix(
    softmax_scale=4, show_progress_bar=False
)
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
k = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(k)

# %%
# Next, we need to compute the Schur vectors and the macrostates. We refer the reader to
# :ref:`sphx_glr_auto_examples_estimators_compute_macrostates.py` where the method is explained more in detail.
g.compute_schur(n_components=4)
g.compute_macrostates(cluster_key="clusters")

# %%
# For :class:`cellrank.estimators.GPCCA`, there are 3 methods for choosing the terminal states:
#
#     1. :meth:`cellrank.estimators.GPCCA.set_terminal_states`
#     2. :meth:`cellrank.estimators.GPCCA.set_terminal_states_from_macrostates`
#     3. :meth:`cellrank.estimators.GPCCA.compute_terminal_states`
#
# We will cover each of these methods below. In the last 2 cases, parameter ``n_cells`` controls how many cells to take
# from each terminal state we take as a categorical annotation.

# %%
# Set terminal states
# ^^^^^^^^^^^^^^^^^^^
# :meth:`cellrank.estimators.GPCCA.set_terminal_states` simply sets the terminal states manually - this
# can be useful when the terminal states are known beforehand. In this case, we don't need to compute the macrostates.
#
# The states can be specified either as a categorical :class:`pandas.Series` where `NaN` values mark cells
# not belonging to any terminal state or a :class:`dict`, where keys correspond to the names of the terminal states,
# and the values to the sequence of cell names or their indices.
#
# Below we set the terminal state called `"Alpha"` as all the cells from the `"Alpha"`
# cluster under ``adata.obs["clusters"]``.
g.set_terminal_states({"Alpha": adata[adata.obs["clusters"] == "Alpha"].obs_names})

# %%
# Set terminal states from macrostates
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# :meth:`cellrank.estimators.GPCCA.set_terminal_states_from_macrostates` sets the terminal states by subsetting
# the macrostates. Note that multiple states can also be combined into new, joint states, as shown below,
# where we combine `"Alpha"` and `"Beta"` states into a new one.
g.set_terminal_states_from_macrostates(["Alpha, Beta", "Epsilon"])

# %%
# Compute terminal states
# ^^^^^^^^^^^^^^^^^^^^^^^
# Lastly, :meth:`cellrank.estimators.GPCCA.compute_terminal_states` which also makes use of the coarse-grained
# transition matrix :attr:`cellrank.estimators.GPCCA.coarse_T` of the macrostates or the `eigengap`
# statistic.
#
# In the example below, we use ``method='eigenap'`` which selects the number of states based on the `eigengap`. The
# terminal states are defined as the top most likely states from the diagonal of the coarse-grained transition matrix.
# To find out more, see :ref:`sphx_glr_auto_examples_estimators_compute_coarse_T.py`.
g.compute_terminal_states(method="eigengap")

# %%
# Now that the terminal states have been either set or computed, we can visualize them in an embedding.
# All of the options seen in :ref:`sphx_glr_auto_examples_estimators_compute_macrostates.py` also apply here,
# like plotting in the same plot (parameter ``same_plot``) or plotting the discrete values (parameter
# ``discrete``).
g.plot_terminal_states(same_plot=False)

# %%
# We note that membership degree of macrostates/terminal states should not be confused with the probability of
# traveling/developing towards these states. For that, we compute the absorption probabilities, see
# :ref:`sphx_glr_auto_examples_estimators_compute_abs_probs.py`. The assignment of cells to macrostates is
# a soft assignment that specifies the degree of membership of any particular cell to a given state.
