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
k = cr.tl.transition_matrix(adata, show_progress_bar=False)
g = cr.tl.estimators.GPCCA(k)

# %%
# First, we need to compute the Schur vectors. By default, only the first 10 vectors are computed.
g.compute_schur()

# %%
# Now we can compute the metastable states of the Markov chain. By default, the number of states is estimated
# by the `eigengap`, which is computed when computing the Schur vectors. Important parameter is the ``cluster_key``,
# which tries to associate the names of metastable states with cluster labels saved under that key.
#
# The second import parameter is ``n_cells``, which select the top number of cells from each state. We are using the
# default 30 cells here.
g.compute_metastable_states(cluster_key="clusters")

# %%
# After computing the metastable states, we can easily inspect them as follows. Below we show the probability of each
# cell belonging to the computed metastable states.
g.metastable_states_probabilities

# %%
# And to get the most probable cells for each state, we can inspect the attribute below.
g.metastable_states

# %%
# Other option how to specify the number of states is to use the `minChi` criterion from from [GPCCA18]_.
# In order to do this, we need to supply a closed interval where we expect our number of metastable states to lie.
g.compute_metastable_states(n_states=[3, 6], use_min_chi=True, cluster_key="clusters")
g.metastable_states_probabilities

# %%
# In the case above, the number of states has stayed the same.
#
# Finally, we can plot the probabilities, as well as the categorical assignment.
g.plot_metastable_states()
g.plot_metastable_states(discrete=True)

# %%
# Both of these options are shown in the same plot, which is not necessarily. To change this, simply run the following.
g.plot_metastable_states(same_plot=False)

# %%
# Lastly, it's also possible to look only at a subset of metastable states or to combine these states into a new one.
# Below we combine the `'Alpha'` and  `'Beta'` - each cell's value corresponds to the probability of belonging
# to the `'Alpha'` or the `'Beta'` metastable state.
g.plot_metastable_states(lineages=["Alpha, Beta"])
