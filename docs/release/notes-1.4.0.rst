CellRank 1.4.0 (2021-06-30)
===========================

Bugfixes
--------

- Fix estimator's inconsistent state when reading from :class:`anndata.AnnData`
  `PR 563 <https://github.com/theislab/cellrank/pull/563>`_.
- Fix not checking whether probabilities sum to 1 in
  :meth:`cellrank.tl.estimators.BaseEstimator.compute_absorption_probabilities`
  `PR 566 <https://github.com/theislab/cellrank/pull/566>`_.
- Fix always forcing sparse transition matrix in :class:`cellrank.tl.kernels.Kernel`
  `PR 586 <https://github.com/theislab/cellrank/pull/586>`_.
- Fix passing custom connectivity key in :class:`cellrank.tl.kernels.Kernel`
  `PR 590 <https://github.com/theislab/cellrank/pull/590>`_.
- Fix kernels in :mod:`cellrank.external` always requiring connectivities
  `PR 600 <https://github.com/theislab/cellrank/pull/600>`_.
- Fix parallelization of sparse matrix with too many jobs `PR 633 <https://github.com/theislab/cellrank/pull/633>`_.
- Fix plotting coarse-grained transition matrix when no stationary distribution is found
  `Issue 594 <https://github.com/theislab/cellrank/issues/594>`_.

Features
--------

- Add new external kernel based on Waddington optimal transport :cite:`schiebinger:19`
  :class:`cellrank.external.kernels.WOTKernel` `PR 595 <https://github.com/theislab/cellrank/pull/595>`_.
- Add new tutorial `Time series datasets <https://cellrank.readthedocs.io/en/stable/real_time.html>`_ that
  shows how to use the external Waddington OT kernel `PR 35 <https://github.com/theislab/cellrank_notebooks/pull/35>`_.
- Add reprogramming dataset :func:`cellrank.datasets.reprogramming_schiebinger` from :cite:`schiebinger:19`
  `PR 631 <https://github.com/theislab/cellrank/pull/631>`_.
- Add :func:`cellrank.pl.log_odds` to plot log-odds ratio between 2 lineages sorted by experimental time.
  `PR 642 <https://github.com/theislab/cellrank/pull/642>`_.
- Add :meth:`cellrank.tl.estimators.GPCCA.plot_macrostate_composition` to visualize macrostate composition over a
  categorical variable `PR 641 <https://github.com/theislab/cellrank/pull/641>`_.
- Add :meth:`cellrank.tl.estimators.BaseEstimator.plot_lineage_drivers_correlation` to plot lineage driver correlation
  with 2 lineages in a scatter plot `PR 640 <https://github.com/theislab/cellrank/pull/640>`_.
- Add :meth:`cellrank.tl.kernels.Kernel.plot_single_flow` based on :cite:`mittnenzweig:21` to visualize outgoing
  transition matrix flow in clusters across experimental time `PR 615 <https://github.com/theislab/cellrank/pull/615>`_.
- Dramatically speed-up the computation of fate probabilities (>5x) `PR 638 <https://github.com/theislab/cellrank/pull/638>`_.
- Remove 4 technical examples and add 2 new examples `PR 602 <https://github.com/theislab/cellrank/pull/602>`_:

  - projection of transition matrix onto an embedding.
  - simulation of random walks on a Markov chain.

- Add option to visualize cell-level covariates in :func:`cellrank.pl.cluster_lineage`
  `PR 634 <https://github.com/theislab/cellrank/pull/634>`_.
- Add option to force-recompute transition matrix in :func:`cellrank.tl.initial_states` and
  :func:`cellrank.tl.terminal_states` `PR 577 <https://github.com/theislab/cellrank/pull/577>`_.
- Change :class:`cellrank.tl.kernels.PseudotimeKernel` defaults and prune available parameters
  of soft threshold scheme `PR 583 <https://github.com/theislab/cellrank/pull/583>`_.
- Parallelize transition matrix computation in :class:`cellrank.tl.kernels.PseudotimeKernel`
  `PR 587 <https://github.com/theislab/cellrank/pull/587>`_.
- Prune *requirements.txt* `PR 571 <https://github.com/theislab/cellrank/pull/571>`_.
- Add small improvements to documentation `PR 584 <https://github.com/theislab/cellrank/pull/584>`_
  `PR 601 <https://github.com/theislab/cellrank/pull/601>`_ `PR 605 <https://github.com/theislab/cellrank/pull/605>`_
  `PR 648 <https://github.com/theislab/cellrank/issues/648>`_.
