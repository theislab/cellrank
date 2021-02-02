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
    :members:

Correlation scheme
++++++++++++++++++

.. autoclass:: cellrank.tl.kernels.CorrelationScheme
    :members:

Dot product scheme
++++++++++++++++++

.. autoclass:: cellrank.tl.kernels.DotProductScheme
    :members:


Connectivity Kernel
-------------------

.. autoclass:: cellrank.tl.kernels.ConnectivityKernel
    :members:

Palantir Kernel
---------------

.. autoclass:: cellrank.tl.kernels.PalantirKernel
    :members:

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


BaseModel
---------

.. autoclass:: cellrank.ul.models.BaseModel
    :members:

Lineage
-------

.. autoclass:: cellrank.tl.Lineage
    :members: reduce, plot_pie, entropy, X, T, view, names, colors
