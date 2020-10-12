# -*- coding: utf-8 -*-
"""
Fit models and plot gene trends
-------------------------------

This example shows how to prepare, fit and plot various models to gene expression trends.

We will focus mostly on Generalized Additive Models (`GAMs <https://en.wikipedia.org/wiki/Generalized_additive_model>`_)
and show how to do this for :mod:`sklearn` estimators towards the end. GAMs are flexible models that are well suited to
model non-linear gene trends as they often appear in single-cell data. Further, they have the advantage that it is
relatively straightforward to derive a confidence interval around the main trend.
"""

from sklearn.svm import SVR

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we need to compute a set of terminal states and we have to estimate fate probabilities towards these. The fate
# probabilities will be used as weights when fitting the model - they determine how much each cell contributes to each
# lineage.
cr.tl.terminal_states(
    adata,
    cluster_key="clusters",
    weight_connectivities=0.2,
    softmax_scale=4,
    show_progress_bar=False,
)
cr.tl.lineages(adata)

# %%
# Models in :mod:`cellrank.ul.models` follow similar patterns as :mod:`sklearn` models.
# We begin by initializing and preparing the model for fitting. :meth:`cellrank.ul.models.BaseModel.prepare` requires
# only the gene name and the lineage name and must be called before :meth:`cellrank.ul.models.BaseModel.fit`. It also
# includes various useful parameters, such as ``time_range`` or ``weight_threshold``, which determine the start and
# end pseudotime and the minimum required threshold for lineage probabilities, respectively.
model = cr.ul.models.GAM(adata)
model.prepare(
    gene="Pak3",
    lineage="Alpha",
    time_key="dpt_pseudotime",
    data_key="Ms",
    n_test_points=100,
)
# %%
# CellRank allows any pseudotime from ``adata.obs`` to be passed via the ``time_key``, e.g. Diffusion Pseudotime (DPT)
# [Haghverdi16]_, :mod:`scvelo`’s latent time [Bergen20]_, Palantir’s pseudotime [Setty19]_, etc.
#
# Further, CellRank accepts imputed gene expression values stored in ``adata.layers`` via the ``data_key``, i.e. you
# can pass MAGIC [MAGIC18]_ imputed data, :func:`scvelo.pp.moments` [Bergen20]_ (used below) or any other form of
# imputation.

# %%
# Once the model has been prepared, it is ready for fitting and prediction.
y_hat = model.fit().predict()
y_hat

# %%
# Optionally, we can also get the confidence interval. Models which don't have a method to compute it, such as
# :class:`cellrank.ul.models.SKLearnModel` wrapper for some :mod:`sklearn` estimators, can use the default the
# confidence interval :meth:`cellrank.ul.models.BaseModel.default_confidence_interval`.
conf_int = model.confidence_interval()
conf_int[:5]

# %%
# After the prediction and optionally the confidence interval calculation, we can plot the smoothed gene expression.
#
# Cells in this plot have been colored by their fate probability of reaching the `Alpha` terminal state. We
# include these probabilities as weights in the loss function when fitting the model. This allows us to weight each
# cell by its relative contribution to the lineage, without needing to subset cells for each lineage.
model.plot(conf_int=True)

# %%
# Lastly, wrapping :mod:`sklearn` estimators is fairly simple, we just pass the instance
# to :class:`cellrank.ul.models.SKLearnModel`.
svr = SVR()
model = cr.ul.models.SKLearnModel(adata, model=svr)
model
