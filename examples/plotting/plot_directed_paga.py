# -*- coding: utf-8 -*-
"""
Plot directed paga
------------------

This example shows how to compute and plot :mod:`scvelo`'s directed PAGA.
"""

import cellrank as cr
import scvelo as scv

adata = cr.datasets.pancreas_preprocessed()

# %%
# First, we compute the initial and terminal states probabilities.
cr.tl.terminal_states(adata, cluster_key="clusters", weight_connectivities=0.2)
cr.tl.initial_states(adata, cluster_key="clusters")

# %%
# We can use :func:`scvelo.tl.recover_latent_time` to compute gene-shared latent time, as well as using initial and
# terminal states computed above. This will be used as a time prior in the computation of the directed PAGA graph.
scv.tl.recover_latent_time(
    adata, root_key="initial_states_probs", end_key="terminal_states_probs"
)


# %%
# Afterwards, we compute the directed PAGA, again specifying te initial and terminal states, as well as the time
# prior mentioned above.
scv.tl.paga(
    adata,
    groups="clusters",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="velocity_pseudotime",
)

# %%
# Lastly, we can plot the results using :func:`cellrank.pl.cluster_fates`.
cr.pl.cluster_fates(
    adata,
    mode="paga_pie",
    cluster_key="clusters",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA",
)
