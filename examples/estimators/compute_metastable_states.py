# -*- coding: utf-8 -*-
"""
Compute metastable states
-------------------------

This examples show how to compute and plot metastable states obtained by :class:`cellrank.tl.estimators.GPCCA`.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, let us prepare the kernel using high-level pipeline and the :class:`cellrank.tl.estimators.GPCCA` estimator.
k = cr.tl.transition_matrix(
    adata, weight_connectivities=0.2, softmax_scale=4, show_progress_bar=False
)
g = cr.tl.estimators.GPCCA(k)

# %%
# First, we need to compute the Schur vectors. By default, only the first 10 vectors are computed.
g.compute_schur()

# %%
# Now we can compute the metastable states of the Markov chain. By default, the number of states is estimated
# by the `eigengap`, which is computed when computing the Schur vectors. Important parameter is the ``cluster_key``,
# which tries to associate the names of metastable states with the cluster labels.
#
# The second import parameter is ``n_cells``, which selects the top cells from each state based
# on the membership degree. By default, only 30 cells are selected.
g.compute_metastable_states(cluster_key="clusters")

# %%
# After computing the metastable states, we can easily inspect them as follows. Below we show for each cell
# the membership degree with respect to the metastable states.
g.metastable_states_probabilities

# %%
# To get the categorical observations, top ``n_cells`` for each metastable state, we can inspect the attribute below.
g.metastable_states

# %%
# Other option how to specify the number of states is to use the `minChi` criterion from [GPCCA18]_.
# In order to do this, we need to supply a closed interval where we expect our number of metastable states to lie.
g.compute_metastable_states(n_states=[3, 6], use_min_chi=True, cluster_key="clusters")
g.metastable_states_probabilities

# %%
# In the case above, the number of states has stayed the same.
#
# Finally, we can plot the membership degrees, as well as the categorical assignment.
g.plot_metastable_states()
g.plot_metastable_states(discrete=True)

# %%
# Both of these options are shown in the same plot, which is not what we always want.
# To change this, simply run the following.
g.plot_metastable_states(same_plot=False)

# %%
# Lastly, it's also possible to look only at a subset of metastable states or to combine these states into new ones.
# Below we combine the `'Alpha'` and  `'Beta'` into joint `'Alpha'` or `'Beta'` metastable state.
g.plot_metastable_states(lineages=["Alpha, Beta"])

# %%
# Method :meth:`cellrank.tl.estimators.GPCCA.compute_metastable_states` also computes the coarse-grained transition
# matrix between metastable states, see :ref:`sphx_glr_auto_examples_estimators_compute_coarse_T.py`.
