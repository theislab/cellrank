API
===

Import CellRank as::

    import cellrank as cr

Once velocities and the velocity graph have been computed using either `scvelo`_ or `velocyto`_,
CellRank offers two modes to interact with its core functionality:

 - high level mode, essentially calling :func:`cellrank.tl.final_states`, :func:`cellrank.tl.root_states` and :func:`cellrank.tl.lineages`. See our `high level tutorial <https://cellrank.readthedocs.io/en/latest/pancreas_basic.html>`_
 - low level mode, interacting directly with the kernels defined in :class:`cellrank.tl.kernels.Kernel` and the estimators :class:`cellrank.tl.GPCCA` or :class:`cellrank.tl.CFLARE`. The division into kernels and estimators ensures that CellRank in broadly applicable, no matter how you have computed your transition matrix.  See our `low level tutorial <https://cellrank.readthedocs.io/en/latest/pancreas_advanced.html>`_.

Additionally, there is a set of plotting functions which can be used downstream of either analysis mode.

The utilities are mainly for fitting continuous models to gene expression data
and are utilised in some of the plotting functions, like :func:`cellrank.pl.gene_trends`.

Tools
~~~~~

.. module:: cellrank.tl
.. currentmodule:: cellrank

.. autosummary::

    :toctree: api/tl

    tl.GPCCA
    tl.CFLARE
    tl.partition
    tl.root_states
    tl.final_states
    tl.lineages
    tl.gene_importance
    tl.transition_matrix
    tl.kernels.VelocityKernel
    tl.kernels.ConnectivityKernel
    tl.kernels.PalantirKernel
    tl.kernels.PrecomputedKernel

Plotting
~~~~~~~~

.. module:: cellrank.pl
.. currentmodule:: cellrank

.. autosummary::

    :toctree: api/pl

    pl.gene_trends
    pl.heatmap
    pl.cluster_lineage
    pl.cluster_fates
    pl.similarity_plot
    pl.graph
    pl.lineages
    pl.composition

Utilities
~~~~~~~~~

.. module:: cellrank.ul
.. currentmodule:: cellrank

.. autosummary::

    :toctree: api/ul

    ul.models.GAM
    ul.models.SKLearnModel
    ul.models.GAMR

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

    :toctree: api/datasets

    datasets.pancreas


.. _scvelo: https://scvelo.readthedocs.io/
.. _velocyto: http://velocyto.org/
