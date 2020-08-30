# -*- coding: utf-8 -*-
"""
Composition plot
================

This plot is used to visualize the composition of a categorical variable within an :class:`anndata.AnnData` object.
"""

import cellrank as cr

# %%
# Load data
# ---------
adata = cr.datasets.pancreas_preprocessed()
adata

# %%
# Plot cluster composition
# ------------------------
cr.pl.composition(adata, key="clusters")
