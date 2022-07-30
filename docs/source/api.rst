API
===
Import CellRank as::

    import cellrank as cr

Once velocities and the velocity graph have been computed using either `scvelo`_ or `velocyto`_,
CellRank offers two modes to interact with its core functionality:

- interacting directly with the kernels defined in :class:`cellrank.kernels.Kernel` and the
  estimators :class:`cellrank.estimators.GPCCA` or :class:`cellrank.estimators.CFLARE`.
  The division into kernels and estimators ensures that CellRank in broadly applicable, no matter how you have
  computed your transition matrix.
  See our `Kernels and estimators tutorial <https://cellrank.readthedocs.io/en/stable/kernels_and_estimators.html>`_.

Additionally, there is a set of plotting functions which can be used downstream of either analysis mode.

The utilities are mainly for fitting continuous models to gene expression data
and are utilized in some of the plotting functions, like :func:`cellrank.pl.gene_trends`.

.. _kernels:

Kernels
-------
Kernels are part of the low-level API and are used to estimate cell-to-cell transitions.

.. module:: cellrank.kernels
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    kernels.VelocityKernel
    kernels.ConnectivityKernel
    kernels.PseudotimeKernel
    kernels.CytoTRACEKernel
    kernels.PrecomputedKernel

External Kernels
----------------
.. module:: cellrank.external
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    external.kernels.StationaryOTKernel
    external.kernels.WOTKernel

.. _estimators:

Estimators
----------
Estimators predict cell fates using the transitions derived from :ref:`Kernels`.

.. module:: cellrank.estimators
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    estimators.GPCCA
    estimators.CFLARE

Plotting
--------
.. module:: cellrank.pl
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    pl.circular_projection
    pl.gene_trends
    pl.log_odds
    pl.heatmap
    pl.cluster_trends
    pl.aggregate_absorption_probabilities

Models
------
.. module:: cellrank.models
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    models.GAM
    models.GAMR
    models.SKLearnModel

Datasets
--------
.. module:: cellrank.datasets
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    datasets.pancreas
    datasets.lung
    datasets.reprogramming_morris
    datasets.reprogramming_schiebinger
    datasets.zebrafish
    datasets.pancreas_preprocessed
    datasets.bone_marrow

.. _scvelo: https://scvelo.readthedocs.io/
.. _velocyto: http://velocyto.org/
