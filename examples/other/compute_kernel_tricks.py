# -*- coding: utf-8 -*-
"""
Kernel tricks
-------------

This example shows some niche, but useful functionalities of :class:`cellrank.tl.kernels.Kernel`.

CellRank is split into :mod:`cellrank.tl.kernels` and :mod:`cellrank.tl.estimators:. Kernels compute transition matrices
based on some inputs, like RNA velocity [Bergen20]_, [Manno18]_, while estimators perform inference based on a given
kernel, e.g. they compute initial and terminal cells and fate probabilities.

Here, will will dive a bit deeper into the how these kernel objects work.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %%
# First, we create some kernels which will be used to compute the cell-to-cell transition matrix:
#
# - :class:`cellrank.tl.kernels.ConnectivityKernel` computes the transition matrix using the KNN graph from
#   :func:`scanpy.pp.neighbors` [Wolf18]_. Note that this kernel is by itself directionless and
#   should be used in conjunction with e.g. :class:`cellrank.tl.kernels.VelocityKernel`.
# - :class:`cellrank.tl.kernels.VelocityKernel` uses RNA velocity for the transition matrix computation [Manno18]_
#   [Bergen20]_, but can also take into account uncertainty in RNA velocity.
# - :class:`cellrank.tl.kernels.PalantirKernel` works similarly as in Palantir [Setty19]_ - it orients the edges of
#   the KNN graph constructed in the expression space using the pseudotemporal ordering of cells, such as
#   Diffusion Pseudotime (DPT) [Haghverdi16]_.
ck = cr.tl.kernels.ConnectivityKernel(adata)
vk = cr.tl.kernels.VelocityKernel(adata)
pk = cr.tl.kernels.PalantirKernel(adata)

# %%
# Kernels can be easily combined with each other. All of the constants within parentheses
# (if they are needed) will be normalized to 1 automatically. This effect will be shown after the transition matrix
# has been computed.
10 * vk + 15 * ck

# %%
# We can also build more complex kernel expression. Note that multiple usage of the same kernel instance in the
# expression caches it's transition matrix, so it's only computed once.
k = ((vk + vk + 42 * vk) + ck * 2) + ck * (3 * vk + 4 * pk)
k

# %%
# For complex expression like the one above, one can get the unique base kernels as follows.
k.kernels

# %%
# Kernels can also be inverted (i.e. the ``backward`` attribute is set to the opposite value to the current one).
# Note that this operation does *NOT* create a copy and modifies the kernel inplace, most importantly, the computed
# transition matrix is removed. This makes it more easier to recompute the transition matrix of kernels,
# since the handles point to the original objects.
#
# Note that in the 2nd :func:`print` statement, we access the private attribute - that's because accessing
# :paramref:`cellrank.tl.kernels.Kernel.transition_matrix` computes the transition matrix with default values.
# This happens only with basic kernels and not the kernel expressions and only if they are not part of a
# larger expression.
ck.compute_transition_matrix()
print(ck.transition_matrix is not None)

inv_ck = ~ck
print(ck._transition_matrix is not None)
print(inv_ck is ck)

# %%
# We can also easily copy the kernel objects.
vk.compute_transition_matrix(
    mode="deterministic", softmax_scale=4, show_progress_bar=False
)
print(vk.transition_matrix is not None)

inv_vk = ~(vk.copy())
print(vk._transition_matrix is not None)
print(vk is inv_vk)

# %%
# After the transition matrix is computed, we can inspect the parameters used in its computation.
ck.compute_transition_matrix()
ck.params

# %%
# The transition matrix can be written to :class:`anndata.AnnData` object. It is saved in the ``.obsp`` attribute.
# This also writes the parameters to the ``.uns`` attribute.
ck.write_to_adata(key="transition_matrix")
adata.obsp["transition_matrix"]

# %%
# Precomputed kernels are useful if you have a transition matrix that was computed outside of CellRank that you would
# like to analyze using one of our estimators. It’s essentially an interface to CellRank’s estimators.
#
# Below we supply a transition matrix saved in ``adata.obsp["transition_matrix"]``, but :class:`anndata.AnnData` object
# is not required and we can just supply either :mod:`numpy` or :mod:`scipy.sparse` array. In that case, a minimal
# empty :class:`anndata.AnnData` object is created, as shown below.
pk = cr.tl.kernels.PrecomputedKernel(adata.obsp["transition_matrix"])
pk.adata
