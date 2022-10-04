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

.. toctree::
    :caption: API
    :maxdepth: 1

    kernels
    estimators
    models
    plotting
    datasets
    developer

.. _scvelo: https://scvelo.readthedocs.io/
.. _velocyto: http://velocyto.org/
