Developer API
=============
Under construction.

Kernels
-------
.. autoclass:: cellrank.kernels.Kernel
    :members:
    :inherited-members:

.. autoclass:: cellrank.kernels.UnidirectionalKernel
    :members:
    :inherited-members:

.. autoclass:: cellrank.kernels.BidirectionalKernel
    :members:
    :inherited-members:

Similarity
~~~~~~~~~~
.. autoclass:: cellrank.kernels.utils.SimilarityABC
    :members:
    :special-members: __call__
    :inherited-members:

.. autoclass:: cellrank.kernels.utils.Cosine
    :members: __call__, hessian

.. autoclass:: cellrank.kernels.utils.Correlation
    :members: __call__, hessian

.. autoclass:: cellrank.kernels.utils.DotProduct
    :members: __call__, hessian

Threshold Scheme
~~~~~~~~~~~~~~~~
.. autoclass:: cellrank.kernels.utils.ThresholdSchemeABC
    :members:
    :special-members: __call__
    :inherited-members:

.. autoclass:: cellrank.kernels.utils.HardThresholdScheme
    :members:
    :special-members: __call__

.. autoclass:: cellrank.kernels.utils.SoftThresholdScheme
    :members:
    :special-members: __call__

.. autoclass:: cellrank.kernels.utils.CustomThresholdScheme
    :members:
    :special-members: __call__

Estimators
----------
.. autoclass:: cellrank.estimators.BaseEstimator
    :members:

.. autoclass:: cellrank.estimators.TermStatesEstimator
    :members:

Models
------
.. autoclass:: cellrank.models.BaseModel
    :members:

Lineage
-------
.. autoclass:: cellrank.Lineage
    :members: priming_degree, reduce, plot_pie, from_adata, X, T, view, names, colors
