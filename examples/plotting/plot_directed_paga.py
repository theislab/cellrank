# -*- coding: utf-8 -*-
"""
Plot directed PAGA
------------------

This example shows how to compute and plot directed PAGA [Wolf19_].
"""

import cellrank as cr
import scvelo as scv

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the initial and terminal states probabilities, as well as the absorption probabilities towards
# the terminal states.
cr.tl.terminal_states(
    adata,
    cluster_key="clusters",
    weight_connectivities=0.2,
    n_states=3,
    softmax_scale=4,
    show_progress_bar=False,
)
cr.tl.lineages(adata)

cr.tl.initial_states(adata, cluster_key="clusters", n_states=1, softmax_scale=4)

# %%
# We can use :func:`scvelo.tl.recover_latent_time` to compute gene-shared latent time leveraging the initial and
# terminal states computed above. This will be used as a time prior when computing the directed PAGA graph.
scv.tl.recover_latent_time(
    adata, root_key="initial_states_probs", end_key="terminal_states_probs"
)

# %%
# Afterwards, we compute the directed PAGA using :func:`scvelo.tl.paga` by again specifying the initial
# and terminal states and the time prior.
scv.tl.paga(
    adata,
    groups="clusters",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="velocity_pseudotime",
)

# %%
# Lastly, we can plot the results using :func:`cellrank.pl.cluster_fates`. For more option, see
# :ref:`sphx_glr_auto_examples_plotting_plot_cluster_fates.py`.
cr.pl.cluster_fates(
    adata,
    mode="paga_pie",
    cluster_key="clusters",
    basis="umap",
    legend_kwargs={"loc": "bottom right"},
    legend_loc="on data",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA",
)
