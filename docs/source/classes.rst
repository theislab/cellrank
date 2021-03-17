Classes
=======

Estimators
~~~~~~~~~~

GPCCA
-----
.. autoclass:: cellrank.tl.estimators.GPCCA
    :members:
    :inherited-members:

CFLARE
------

.. autoclass:: cellrank.tl.estimators.CFLARE
    :members:
    :inherited-members:

Kernels
~~~~~~~

Velocity Kernel
---------------

.. autoclass:: cellrank.tl.kernels.VelocityKernel
    :members:

Cosine similarity scheme
++++++++++++++++++++++++

.. autoclass:: cellrank.tl.kernels.CosineScheme
    :members: __call__, hessian

Correlation scheme
++++++++++++++++++

.. autoclass:: cellrank.tl.kernels.CorrelationScheme
    :members: __call__, hessian

Dot product scheme
++++++++++++++++++

.. autoclass:: cellrank.tl.kernels.DotProductScheme
    :members: __call__, hessian


Connectivity Kernel
-------------------

.. autoclass:: cellrank.tl.kernels.ConnectivityKernel
    :members:

Pseudotime Kernel
-----------------

.. autoclass:: cellrank.tl.kernels.PseudotimeKernel
    :members:

Hard threshold scheme
+++++++++++++++++++++

.. autoclass:: cellrank.tl.kernels.HardThresholdScheme
    :members:
    :special-members: __call__

Soft threshold scheme
+++++++++++++++++++++

.. autoclass:: cellrank.tl.kernels.SoftThresholdScheme
    :members:
    :special-members: __call__


Precomputed Kernel
------------------

.. autoclass:: cellrank.tl.kernels.PrecomputedKernel
    :members:

Models
~~~~~~

GAM
---

.. autoclass:: cellrank.ul.models.GAM
    :members:
    :inherited-members:

SKLearnModel
------------

.. autoclass:: cellrank.ul.models.SKLearnModel
    :members:
    :inherited-members:

GAMR
----

.. autoclass:: cellrank.ul.models.GAMR
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

Similarity scheme
-----------------

.. autoclass:: cellrank.tl.kernels.SimilaritySchemeABC
    :members:
    :special-members: __call__
    :inherited-members:


Threshold scheme
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
    :members: reduce, plot_pie, entropy, X, T, view, names, colors
