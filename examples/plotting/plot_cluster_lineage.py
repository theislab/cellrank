"""
Plot cluster lineage
--------------------

This example shows how to cluster and plot genes in a specific lineage.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the absorption probabilities and select a model that will be used for gene trend smoothing.
cr.tl.terminal_states(
    adata,
    cluster_key="clusters",
    weight_connectivities=0.2,
    n_states=3,
    softmax_scale=4,
    show_progress_bar=False,
)
cr.tl.lineages(adata)

model = cr.ul.models.GAM(adata)

# %%
# Next, we can fit the model for some subset of genes for a specific lineage, as seen below. After the model
# has been fitted, we use it to get the smoothed gene expression for the test points (by default, it is 200 points
# uniformly spaced along the pseudotime). Afterwards, we reduce the dimension using `PCA` and cluster using the
# `louvain` algorithm.
#
# Note that calling this function twice will use the already computed values, unless ``recompute=True`` is specified.
cr.pl.cluster_lineage(
    adata,
    model,
    adata.var_names[:200],
    lineage="Alpha",
    time_key="dpt_pseudotime",
    show_progress_bar=False,
)


# %%
# The clustered genes can be accessed as shown below. In general, ``adata.uns['lineage_..._trend']`` contains
# :class:`anndata` object of shape `(n_test_points, n_genes)`.
adata.uns["lineage_Alpha_trend"].obs["clusters"]
