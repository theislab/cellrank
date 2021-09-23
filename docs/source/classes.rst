Classes
=======

Estimators
~~~~~~~~~~

GPCCA
-----
.. autoclass:: cellrank.tl.estimators.GPCCA
    :noindex:
    :members:
    :inherited-members:

CFLARE
------
.. autoclass:: cellrank.tl.estimators.CFLARE
    :noindex:
    :members:
    :inherited-members:

Kernels
~~~~~~~

Velocity Kernel
---------------
.. autoclass:: cellrank.tl.kernels.VelocityKernel
    :noindex:
    :members:

Cosine Similarity Scheme
++++++++++++++++++++++++
.. autoclass:: cellrank.tl.kernels.CosineScheme
    :members: __call__, hessian

Correlation Scheme
++++++++++++++++++
.. autoclass:: cellrank.tl.kernels.CorrelationScheme
    :members: __call__, hessian

Dot Product Scheme
++++++++++++++++++
.. autoclass:: cellrank.tl.kernels.DotProductScheme
    :members: __call__, hessian

Connectivity Kernel
-------------------
.. autoclass:: cellrank.tl.kernels.ConnectivityKernel
    :noindex:
    :members:

Pseudotime Kernel
-----------------
.. autoclass:: cellrank.tl.kernels.PseudotimeKernel
    :noindex:
    :members:

Hard Threshold Scheme
+++++++++++++++++++++
.. autoclass:: cellrank.tl.kernels.HardThresholdScheme
    :members:
    :special-members: __call__

Soft Threshold Scheme
+++++++++++++++++++++
.. autoclass:: cellrank.tl.kernels.SoftThresholdScheme
    :members:
    :special-members: __call__

CytoTRACE Kernel
----------------
.. autoclass:: cellrank.tl.kernels.CytoTRACEKernel
    :noindex:
    :members: compute_cytotrace, compute_transition_matrix

Precomputed Kernel
------------------
.. autoclass:: cellrank.tl.kernels.PrecomputedKernel
    :noindex:
    :members:

Models
~~~~~~

GAM
---
.. autoclass:: cellrank.ul.models.GAM
    :noindex:
    :members:
    :inherited-members:

SKLearnModel
------------
.. autoclass:: cellrank.ul.models.SKLearnModel
    :noindex:
    :members:
    :inherited-members:

GAMR
----
.. autoclass:: cellrank.ul.models.GAMR
    :noindex:
    :members:
    :inherited-members:

Base Classes
~~~~~~~~~~~~

BaseEstimator
-------------
.. autoclass:: cellrank.tl.estimators.BaseEstimator
    :members:

Kernel
------
.. autoclass:: cellrank.tl.kernels.Kernel
    :members:
    :inherited-members:

ExperimentalTime Kernel
-----------------------
.. autoclass:: cellrank.tl.kernels.ExperimentalTimeKernel
    :members:
    :inherited-members:

TransportMap Kernel
-------------------
.. autoclass:: cellrank.tl.kernels.TransportMapKernel
    :members:
    :inherited-members:

Similarity Scheme
-----------------
.. autoclass:: cellrank.tl.kernels.SimilaritySchemeABC
    :members:
    :special-members: __call__
    :inherited-members:

Threshold Scheme
----------------
.. autoclass:: cellrank.tl.kernels.ThresholdSchemeABC
    :members:
    :special-members: __call__
    :inherited-members:

BaseModel
---------
.. autoclass:: cellrank.ul.models.BaseModel
    :members:

Lineage
-------
.. autoclass:: cellrank.tl.Lineage
    :members: priming_degree, reduce, plot_pie, from_adata, X, T, view, names, colors
