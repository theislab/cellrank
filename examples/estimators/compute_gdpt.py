# -*- coding: utf-8 -*-
"""
Compute gDPT
------------

This examples shows how to compute and plot generalized Diffusion pseudotime from [Haghverdi16]_ .
"""

import cellrank as cr
import scanpy as sc
import scvelo as scv
import numpy as np

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata.uns["iroot"] = root_idx = np.where(adata.obs["clusters"] == "Ngn3 low EP")[0][0]
sc.tl.dpt(adata)
adata

# %%
# First, let us prepare the kernel using high-level pipeline and the :class:`cellrank.tl.estimators.GPCCA` estimator.
k = cr.tl.transition_matrix(adata, show_progress_bar=False)
g = cr.tl.estimators.GPCCA(k)

# %%
# Generalized DPT is exactly the same as the original, except that instead of eigenvectors,
# we use the real Schur vectors.
#
# Note that the line below does not write any attribute to our estimator object, only writes to ``adata.obs``
# based on parameter ``key_added`` (defaults to `'gdpt_pseudotime'`).
g.compute_gdpt(n_components=20)

# %%
# We can now compare the with the original Diffusion pseudotime from :mod:`scanpy` and the generalized version.
scv.pl.scatter(adata, color=["dpt_pseudotime", "gdpt_pseudotime"])
