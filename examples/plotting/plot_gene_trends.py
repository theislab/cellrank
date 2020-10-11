# -*- coding: utf-8 -*-
"""
Plot gene trends
----------------

This example shows how to plot smoothed gene expression toward specific terminal populations.

By default, we use Generalized Additive Models (`GAMs <https://en.wikipedia.org/wiki/Generalized_additive_model>`_)
to fit gene expression values and we specify each cell’s contribution to each lineage via the lineage probabilities.

For models based on :mod:`sklearn` estimators, see :class:`cellrank.ul.models.SKLearnModel`.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we compute the terminal states and the absorption probabilities towards them.
# The absorption probabilities will be used as weights in the loss function when fitting the GAMs.
#
# We further set up a model instance that will be used for smoothing.
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
# To plot the trends for some genes, run the code below. The parameter ``data_key`` specifies layer in ``adata.layers``
# from which we obtain gene expression - in this case, we use moments of gene expression computed by
# :func:`scvelo.pp.moments`.
cr.pl.gene_trends(
    adata,
    model,
    ["Map2", "Dcx"],
    data_key="Ms",
    time_key="dpt_pseudotime",
    show_progress_bar=False,
)

# %%
# We can plot all trends in the same plot and hide the cells, as shown below.
cr.pl.gene_trends(
    adata,
    model,
    ["Map2", "Dcx"],
    data_key="Ms",
    same_plot=True,
    hide_cells=True,
    time_key="dpt_pseudotime",
    show_progress_bar=False,
)

# %%
# We can specify specific pseudotime ranges, separately for each lineage - the model is always fitted
# using cells with the given time range. By default, the start is always 0 and the end is automatically determined
# based on the absorption probabilities.
cr.pl.gene_trends(
    adata,
    model,
    ["Map2", "Dcx"],
    data_key="Ms",
    lineages=["Alpha", "Beta"],
    time_range=[(0.2, 1), (0, 0.8)],
    time_key="dpt_pseudotime",
    show_progress_bar=False,
)

# %%
# We can also return the models, which can be useful to inspect the fitted models more granularly or when the
# fitting has failed - such models will be returned as :class:`cellrank.ul.models.FailedModel`.
#
# Below we show what would happen if a model were to fail for some arbitrary gene and lineage combination.
failed_model = cr.ul.models.FailedModel(model, exc="This is just a dummy example.")

models = cr.pl.gene_trends(
    adata,
    {
        "Map2": {"Alpha": failed_model, "*": model},
        "Dcx": {"Beta": failed_model, "*": model},
    },
    ["Map2", "Dcx"],
    data_key="Ms",
    time_key="dpt_pseudotime",
    show_progress_bar=False,
    return_models=True,
)
models["Map2"]["Alpha"]


# %%
# There are many more options as how to customize the plot or how to pass additional arguments to
# :meth:`cellrank.ul.models.BaseModel.prepare` and we encourage the reader to take a closer look at the documentation of
# :func:`cellrank.pl.gene_trends` and the above mentioned method.
