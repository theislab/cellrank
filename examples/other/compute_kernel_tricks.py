# -*- coding: utf-8 -*-
"""
Kernel tricks
-------------

This example shows some niche, but useful functionalities of :class:`cellrank.tl.kernels.Kernel`.
"""

import cellrank as cr

adata = cr.datasets.pancreas_preprocessed("../example.h5ad")
vk = cr.tl.kernels.VelocityKernel(adata)
ck = cr.tl.kernels.ConnectivityKernel(adata)

# %%
# Kernels can be easily and freely combined with each other. All of the constants within parentheses
# (if they are needed) will be normalized to 1 automatically.
10 * vk + 15 * ck

# %%
# We can build much more complex kernel expression. Note that multiple invocation of the same kernel instance in the
# expression caches it's transition matrix, so it's only computed once.
k = ((vk + vk + 42 * vk) + ck * 2) + ck * vk
k

# %%
# For complex expression like the one above, one can get the basic kernels as follows.
k.kernels

# %%
# Kernels can also be inverted (i.e. the ``backward`` attribute is set to the opposite value to the current one).
# Note that this operation does NOT create a copy and modifies the kernel inplace, most importantly, the computed
# transition matrix is removed. This makes it slightly more easier to recompute the transition matrix of kernels,
# since the objects are the same.
#
# Note that in the 2nd :func:`print` statement, we access the private attribute -
# that's because accessing :paramref:`cellrank.tl.kernels.Kernel.transition_matrix` computes the transition matrix
# with default values. This happens only with basic kernels, not the kernel expressions.
ck.compute_transition_matrix()
print(ck.transition_matrix is not None)

inv_ck = ~ck
print(ck._transition_matrix is not None)
print(inv_ck is ck)

# %%
# Kernels can be copied, which can solve the problem
vk.compute_transition_matrix(mode="deterministic", show_progress_bar=False)
print(vk.transition_matrix is not None)

inv_vk = ~(vk.copy())
print(vk._transition_matrix is not None)
print(vk is inv_vk)

# %%
# After the transition matrix is computed, we can inspect important parameters used in its computation.
ck.compute_transition_matrix()
ck.params

# %%
# The transition matrix can be written to :mod:`anndata` object. It is saved in the ``.obsp`` attribute and the key
# is optional. This also writes the parameters to ``.uns`` attribute of the :mod:`anndata` object.
ck.write_to_adata(key="transition_matrix")
adata.obsp["transition_matrix"]

# %%
# Precomputed kernel doesn't require :mod:`anndata` object - a dummy one will be created if none is supplied.
# This can be useful in conjunction with :meth:`cellrank.tl.estimators.BaseEstimator.set_final_states` to calculate
# absorption probabilities using :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities`.
pk = cr.tl.kernels.PrecomputedKernel(adata.obsp["transition_matrix"])
pk.adata
