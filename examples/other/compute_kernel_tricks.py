# -*- coding: utf-8 -*-
"""
Kernel tricks
-------------

This example shows some niche, but useful functionalities of :class:`cellrank.tl.kernels.Kernel`.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
adata

# %% First, we create some kernels that will be used to compute the cell-to-cell transition matrix.
vk = cr.tl.kernels.VelocityKernel(adata)
ck = cr.tl.kernels.ConnectivityKernel(adata)
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
# Note that in the 2nd :func:`print` statement, we access the private attribute -
# that's because accessing :paramref:`cellrank.tl.kernels.Kernel.transition_matrix` computes the transition matrix
# with default values. This happens only with basic kernels and not the kernel expressions and only if they are not part
# of a larger expression.
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
# Precomputed kernel doesn't require :mod:`anndata` object - a dummy one will be created if none is supplied.
# This can be useful in conjunction with :meth:`cellrank.tl.estimators.BaseEstimator.set_terminal_states` to calculate
# absorption probabilities using :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities`.
pk = cr.tl.kernels.PrecomputedKernel(adata.obsp["transition_matrix"])
pk.adata
