.. role:: small

1.4.0 :small:`SOON`
~~~~~~~~~~~~~~~~~~~
This release includes:

.. rubric:: Additions

- Add new external kernel based on Waddington optimal transport :cite:`schiebinger:19`
  :class:`cellrank.external.kernels.WOTKernel` `PR 595 <https://github.com/theislab/cellrank/pull/595>`_.
- Remove 4 technical examples and add 2 new examples `PR 602 <https://github.com/theislab/cellrank/pull/602>`_:

  - :ref:`sphx_glr_auto_examples_kernels_plot_projection.py` - transition matrix project onto an embedding.
  - :ref:`sphx_glr_auto_examples_kernels_plot_random_walks.py` - simulation of random walks on a Markov chain.

- Add option to force-recompute transition matrix in :func:`cellrank.tl.initial_states` and
  :func:`cellrank.tl.terminal_states` `PR 577 <https://github.com/theislab/cellrank/pull/577>`_.
- Change :class:`cellrank.tl.kernels.PseudotimeKernel` defaults and prune available parameters
  of soft thresholding scheme `PR 583 <https://github.com/theislab/cellrank/pull/583>`_.
- Parallelize transition matrix computation in :class:`cellrank.tl.kernels.PseudotimeKernel`
  `PR 587 <https://github.com/theislab/cellrank/pull/587>`_.
- Prune *requirements.txt* `PR 571 <https://github.com/theislab/cellrank/pull/571>`_.
- Add small improvements to documentation `PR 584 <https://github.com/theislab/cellrank/pull/584>`_
  `PR 601 <https://github.com/theislab/cellrank/pull/601>`_ `PR 605 <https://github.com/theislab/cellrank/pull/605>`_.

.. rubric:: Bugfixes

- Fix estimator's incosistent state when reading from :class:`anndata.AnnData`
  `PR 563 <https://github.com/theislab/cellrank/pull/563>`_.
- Fix not checking whether probabilities sum to 1 in
  :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities`
  `PR 566 <https://github.com/theislab/cellrank/pull/566>`_.
- Fix always forcing sparse transition matrix in :class:`cellrank.tl.kernels.Kernel`
  `PR 586 <https://github.com/theislab/cellrank/pull/586>`_.
- Fix passing custom connectivities key in :class:`cellrank.tl.kernels.Kernel`
  `PR 590 <https://github.com/theislab/cellrank/pull/590>`_.
- Fix kernels in :mod:`cellrank.external` always requiring connectivities
  `PR 600 <https://github.com/theislab/cellrank/pull/600>`_.
