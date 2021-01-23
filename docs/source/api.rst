API
===

Import CellRank as::

    import cellrank as cr

Once velocities and the velocity graph have been computed using either `scvelo`_ or `velocyto`_,
CellRank offers two modes to interact with its core functionality:

- high level mode, essentially calling :func:`cellrank.tl.terminal_states`, :func:`cellrank.tl.initial_states` and
  :func:`cellrank.tl.lineages`.
  See our `high level tutorial  <https://cellrank.readthedocs.io/en/latest/pancreas_basic.html>`_
- low level mode, interacting directly with the kernels defined in :class:`cellrank.tl.kernels.Kernel` and the
  estimators :class:`cellrank.tl.estimators.GPCCA` or :class:`cellrank.tl.estimators.CFLARE`.
  The division into kernels and estimators ensures that CellRank in broadly applicable, no matter how you have
  computed your transition matrix.
  See our `low level tutorial <https://cellrank.readthedocs.io/en/latest/pancreas_advanced.html>`_.

Additionally, there is a set of plotting functions which can be used downstream of either analysis mode.

The utilities are mainly for fitting continuous models to gene expression data
and are utilized in some of the plotting functions, like :func:`cellrank.pl.gene_trends`.

Tools
~~~~~
This module offers a high-level API to compute cell fates and driver genes.

.. module:: cellrank.tl
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    tl.transition_matrix
    tl.initial_states
    tl.terminal_states
    tl.lineages
    tl.lineage_drivers

.. _kernels:

Kernels
-------
Kernels are part of the low-level API and are used to estimate cell-to-cell transitions.

.. module:: cellrank.tl
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    tl.kernels.VelocityKernel
    tl.kernels.ConnectivityKernel
    tl.kernels.PalantirKernel
    tl.kernels.PrecomputedKernel

Estimators
----------
Estimators predict cell fates using the transitions derived from :ref:`Kernels`.

.. module:: cellrank.tl
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    tl.estimators.GPCCA
    tl.estimators.CFLARE

Plotting
~~~~~~~~

.. module:: cellrank.pl
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    pl.initial_states
    pl.terminal_states
    pl.lineages
    pl.lineage_drivers
    pl.circular_projection
    pl.gene_trends
    pl.heatmap
    pl.cluster_lineage
    pl.cluster_fates
    pl.graph

Utilities
~~~~~~~~~

.. module:: cellrank.ul
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    ul.models.GAM
    ul.models.GAMR
    ul.models.SKLearnModel
    ul.models.FittedModel
    ul.models.FailedModel

Reading
~~~~~~~
.. module:: cellrank
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    read

Datasets
~~~~~~~~

.. module:: cellrank.datasets
.. currentmodule:: cellrank

.. autosummary::
    :toctree: api

    datasets.pancreas
    datasets.lung
    datasets.reprogramming
    datasets.pancreas_preprocessed


.. _scvelo: https://scvelo.readthedocs.io/
.. _velocyto: http://velocyto.org/
