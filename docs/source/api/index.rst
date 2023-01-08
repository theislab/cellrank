API
===
Import CellRank as::

    import cellrank as cr

CellRank has a modular API, organized around kernels and estimators:

- :doc:`Kernels <kernels>` compute cell-cell transition matrices using various input data modalities,
  including RNA velocity, any pseudotime, a developmental potential, experimental time points, and more.
- :doc:`Estimators <estimators>` use the cell-cell transition matrix to derive insights about cellular dynamics,
  for example, they compute initial and terminal states, fate probabilities, and driver genes. Our recommended
  estimator is :class:`cellrank.estimators.GPCCA`.
- In addition, there are :doc:`models <models>` for gene trend fitting, :doc:`plotting functions <plotting>`
  for visualization, and :doc:`datasets <datasets>` that help you getting started with CellRank.

If you are new to CellRank, check out our :doc:`about CellRank <../about/index>` page, and take a look at our
:doc:`tutorials <../tutorials>`.

.. toctree::
    :caption: API
    :maxdepth: 1

    kernels
    estimators
    models
    plotting
    datasets
    developer
