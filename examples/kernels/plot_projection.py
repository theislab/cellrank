"""
Plot projection
---------------

This example shows how to project a transition matrix onto a low-dimensional embedding.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata
