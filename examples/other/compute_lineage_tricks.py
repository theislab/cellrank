# -*- coding: utf-8 -*-
"""
Lineage tricks
--------------

This example shows some niche, but useful functionalities of :class:`cellrank.tl.Lineage`.
"""

import cellrank as cr
import numpy as np

np.random.seed(42)

# %%
# The lineage class behaves like a :mod:`numpy` array, for the most part. The key difference is that
# it tries to always preserve it's 2 dimensional shape and that it has :mod:`pandas`-like indexing.
#
# The constructor requires the underlying array and the lineage names, which must be unique. The colors are optional
# and by default they are automatically generated.

lin = cr.tl.Lineage(
    np.abs(np.random.normal(size=(10, 4))), names=["foo", "bar", "baz", "quux"]
)
lin /= lin.sum(1)

# %%
# In some cases, this behavior is not desirable or can have unintended consequences. To access the underlying
# :class:`numpy` array, use the :paramref:`cellrank.tl.Lineage.X` attribute.
lin.X

# %%
# Lineages can also be transposed.
lin.T

# %%
# Indexing into lineage can be done via the names as well.
lin[["foo", "bar"]]

# %%
# Two or more lineage can be combined into by joining the names with `","`. This also automatically
# updates the color based on the combined lineages' colors.
lin[["bar, baz, quux"]]

# %%
# Most of the :mod:`numpy` methods are supported by the :class:`cellrank.tl.Lineage`. One can also calculate the
# entropy, which in [Setty19]_ is defined as the differentiation potential.
lin.entropy(axis=1)

# %%
# When subsetting the lineage and not selecting all of them, they will no longer sum to 1 and cannot be
# interpreted as a probability distribution. We offer a method :paramref:`cellrank.tl.Lineage.reduce` which
# can be used to solve this issue. Below we show only one out of many normalization techniques.
lin.reduce("foo, quux", "baz", normalize_weights="softmax")

# %%
# Lastly, we can plot aggregate information about lineages, such as :func:`numpy.mean` and others.
lin.plot_pie(np.mean, legend_loc="on data")
