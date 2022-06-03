Classes
=======

Estimators
~~~~~~~~~~

GPCCA
-----
.. autoclass:: cellrank.estimators.GPCCA
    :noindex:
    :members:
    :inherited-members:

CFLARE
------
.. autoclass:: cellrank.estimators.CFLARE
    :noindex:
    :members:
    :inherited-members:

Kernels
~~~~~~~

Velocity Kernel
---------------
.. autoclass:: cellrank.kernels.VelocityKernel
    :noindex:
    :members:

Cosine Similarity Scheme
++++++++++++++++++++++++
.. autoclass:: cellrank.kernels.utils.Cosine
    :members: __call__, hessian

Correlation Scheme
++++++++++++++++++
.. autoclass:: cellrank.kernels.utils.Correlation
    :members: __call__, hessian

Dot Product Scheme
++++++++++++++++++
.. autoclass:: cellrank.kernels.utils.DotProduct
    :members: __call__, hessian

Connectivity Kernel
-------------------
.. autoclass:: cellrank.kernels.ConnectivityKernel
    :noindex:
    :members:

Pseudotime Kernel
-----------------
.. autoclass:: cellrank.kernels.PseudotimeKernel
    :noindex:
    :members:

Hard Threshold Scheme
+++++++++++++++++++++
.. autoclass:: cellrank.kernels.utils.HardThresholdScheme
    :members:
    :special-members: __call__

Soft Threshold Scheme
+++++++++++++++++++++
.. autoclass:: cellrank.kernels.utils.SoftThresholdScheme
    :members:
    :special-members: __call__

CytoTRACE Kernel
----------------
.. autoclass:: cellrank.kernels.CytoTRACEKernel
    :noindex:
    :members: compute_cytotrace, compute_transition_matrix

Precomputed Kernel
------------------
.. autoclass:: cellrank.kernels.PrecomputedKernel
    :noindex:
    :members:

Models
~~~~~~

GAM
---
.. autoclass:: cellrank.models.GAM
    :noindex:
    :members:
    :inherited-members:

SKLearnModel
------------
.. autoclass:: cellrank.models.SKLearnModel
    :noindex:
    :members:
    :inherited-members:

GAMR
----
.. autoclass:: cellrank.models.GAMR
    :noindex:
    :members:
    :inherited-members:

Base Classes
~~~~~~~~~~~~

BaseEstimator
-------------
.. autoclass:: cellrank.estimators.BaseEstimator
    :members:

Kernel
------
.. autoclass:: cellrank.kernels.Kernel
    :members:
    :inherited-members:

ExperimentalTime Kernel
-----------------------
.. autoclass:: cellrank.kernels.ExperimentalTimeKernel
    :members:
    :inherited-members:

TransportMap Kernel
-------------------
.. autoclass:: cellrank.kernels.TransportMapKernel
    :members:
    :inherited-members:

Similarity Scheme
-----------------
.. autoclass:: cellrank.kernels.utils.SimilarityABC
    :members:
    :special-members: __call__
    :inherited-members:

Threshold Scheme
----------------
.. autoclass:: cellrank.kernels.utils.ThresholdSchemeABC
    :members:
    :special-members: __call__
    :inherited-members:

BaseModel
---------
.. autoclass:: cellrank.models.BaseModel
    :members:

Lineage
-------
.. autoclass:: cellrank.Lineage
    :members: priming_degree, reduce, plot_pie, from_adata, X, T, view, names, colors
