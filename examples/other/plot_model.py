# -*- coding: utf-8 -*-
"""
Fit model and plot gene trend
-----------------------------

This example shows how to prepare and fit :class:`cellrank.ul.models.BaseModel` and how to plot the estimated trend.
"""

from sklearn.svm import SVR

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we need to compute the lineages of the forward process.
cr.tl.terminal_states(
    adata,
    cluster_key="clusters",
    weight_connectivities=0.2,
    softmax_scale=4,
    show_progress_bar=False,
)
cr.tl.lineages(adata)

# %%
# Models in :mod:`cellrank.ul.models` follow similar patterns as :mod:`sklearn` models. Below we show an example usage.
# We begin by initializing and preparing the model for fitting. :meth:`cellrank.ul.models.BaseModel.prepare` requires
# only a gene name and a lineage name and must be called before :meth:`cellrank.ul.models.BaseModel.fit`. It also
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
# After the model has been prepared, it is ready for fitting and prediction.
y_hat = model.fit().predict()
y_hat

# %%
# Optionally, we can also get the confidence interval. Models which don't have a method to compute it, such as
# :class:`cellrank.ul.models.SKLearnModel` wrapper for :mod:`sklearn` estimators, use the default confidence interval
# define here :meth:`cellrank.ul.models.BaseModel.default_confidence_interval`.
conf_int = model.confidence_interval()
conf_int[:5]

# %%
# After the prediction and confidence interval calculation, we can plot the results.
model.plot(show_conf_int=True)

# %%
# Wrapping :mod:`sklearn` estimators is simple, just pass the instance to :class:`cellrank.ul.models.SKLearnModel`.
svr = SVR()
model = cr.ul.models.SKLearnModel(adata, model=svr)
model
