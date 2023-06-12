API
===
Import CellRank as::

    import cellrank as cr

CellRank has a modular API, organized around kernels and estimators:

- :mod:`cellrank.kernels` compute cell-cell transition matrices using various input data modalities,
  including RNA velocity, any pseudotime, a developmental potential, experimental time points, and more.
- :mod:`cellrank.estimators` use the cell-cell transition matrix to derive insights about cellular dynamics,
  for example, they compute initial and terminal states, fate probabilities, and driver genes. Our recommended
  estimator is the :class:`~cellrank.estimators.GPCCA` estimator.
- In addition, there are :mod:`cellrank.models` for gene trend fitting, :mod:`plotting <cellrank.pl>`
  for visualization, and :mod:`cellrank.datasets` that help you getting started with CellRank.

If you are new to CellRank, check out our :doc:`about CellRank <../about/index>` page, and take a look at our
:doc:`tutorials <../notebooks/tutorials/index>`.

.. toctree::
    :caption: API
    :maxdepth: 1

    kernels
    estimators
    models
    plotting
    datasets
    developer
