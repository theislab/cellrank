Release Notes
=============

.. role:: small

Version 1.0
-----------

1.3.0 :small:`2021-03-29`
~~~~~~~~~~~~~~~~~~~~~~~~~
This release includes some major additions which make CellRank more applicable with and without RNA velocity
information. In particular, it includes:

.. rubric:: Additions

- Add new kernel :class:`cellrank.tl.kernels.CytoTRACEKernel` which computes cell-cell transition probabilities based
  on the CytoTRACE score :cite:`gulati:20`, a measure of differentiation potential,
  `PR 527 <https://github.com/theislab/cellrank/pull/527>`_.
- Add external API :mod:`cellrank.external` with a stationary optimal transport kernel
  :class:`cellrank.external.kernels.OTKernel` contributed from :cite:`zhang:21`, as well as a
  `contributing guide <https://github.com/theislab/cellrank/blob/master/CONTRIBUTING.rst>`_,
  `PR 522 <https://github.com/theislab/cellrank/pull/522>`_.
- Rename ``cellrank.tl.kernels.PalantirKernel`` to :class:`cellrank.tl.kernels.PseudotimeKernel` and add
  hard thresholding scheme inspired by :cite:`setty:19`, a soft thresholding scheme inspired by :cite:`stassen:21` and
  a custom scheme when computing the transition matrix, see e.g. :class:`cellrank.tl.kernels.SoftThresholdScheme`
  `PR 514 <https://github.com/theislab/cellrank/pull/514>`_.
- Add more flexibility to :class:`cellrank.tl.kernels.ConnectivityKernel`, allowing it to use any cell-cell similarities
  from :attr:`anndata.AnnData.obsp`, such as spatial similarities from :mod:`squidpy` :cite:`palla:21`
  `PR 501 <https://github.com/theislab/cellrank/pull/501>`_.
- Revamp `Pancreas Advanced <https://cellrank.readthedocs.io/en/latest/pancreas_advanced.html>`_ tutorial
  to showcase CellRank's modular structure of kernels and estimators.
  `PR 32 <https://github.com/theislab/cellrank_notebooks/pull/32>`_.
- Add 2 new tutorials:

  - `Beyond RNA velocity <https://cellrank.readthedocs.io/en/latest/beyond_rna_velocity.html>`_: shows how to use
    CellRank when no RNA velocity information is available.
    `PR 32 <https://github.com/theislab/cellrank_notebooks/pull/32>`_
  - `Creating a new kernel <https://cellrank.readthedocs.io/en/latest/creating_new_kernel.html>`_: explains how to
    create your own custom kernel class that estimates cell-cell transition probabilities
    `PR 31 <https://github.com/theislab/cellrank_notebooks/pull/31>`_.

- Add projection of transition matrix onto an embedding :meth:`cellrank.tl.kernels.Kernel.compute_projection`
- Add random walk simulation and visualization in an embedding :meth:`cellrank.tl.kernels.Kernel.plot_random_walks`
  `PR 537 <https://github.com/theislab/cellrank/pull/537>`_.
- Add :meth:`cellrank.tl.Lineage.priming_degree` `PR 502 <https://github.com/theislab/cellrank/pull/502>`_
  which estimates a cell's plasticity/differentiation potential based on ideas by :cite:`setty:19`
  and :cite:`velten:17`.
- Add checks for transition matrix irreducibility `PR 516 <https://github.com/theislab/cellrank/pull/516>`_.
- Add Zebrafish development dataset from :cite:`farrel:18` `PR 539 <https://github.com/theislab/cellrank/pull/539>`_.
- Speed-up stationary distribution calculation in :mod:`pygpcca` `PR 22 <https://github.com/msmdev/pyGPCCA/pull/22>`_.

.. rubric:: Bugfixes

- Fix various bugs when plotting multiple gene trends `PR 487 <https://github.com/theislab/cellrank/pull/487>`_.
- Fix gene trend smoothing not working for 1 lineage `PR 512 <https://github.com/theislab/cellrank/pull/512>`_.
- Fix :mod:`pandas` error when computing macrostates `PR 513 <https://github.com/theislab/cellrank/pull/513>`_.
- Remove malfunctioning *Edit on GitHub* from the documentation
  `PR 538 <https://github.com/theislab/cellrank/pull/538>`_.

1.2.0 :small:`2021-02-02`
~~~~~~~~~~~~~~~~~~~~~~~~~
This release includes:

.. rubric:: Additions

- Completely **refactored the underlying code base of GPCCA** and set it up as it's own package called
  `pyGPCCA <https://pygpcca.readthedocs.io/en/latest/>`_ with documentation and an example. Going forwards, this will
  ensure that one of the "engines" of CellRank is also easy to maintain to extend. Further, this will make CellRank's
  installation more convenient by not needing to vendorize additional dependencies
  `PR 472 <https://github.com/theislab/cellrank/pull/472>`_.
- Add :func:`cellrank.pl.circular_projection` visualizing computed fate probabilities as done in :cite:`velten:17`,
  see :ref:`sphx_glr_auto_examples_plotting_plot_circular_embedding.py`.
  `PR 459 <https://github.com/theislab/cellrank/pull/459>`_.
- Allow legends not to be plotted by passing ``legend_loc="none"``, as done in `scVelo <https://scvelo.org>`_
  `PR 470 <https://github.com/theislab/cellrank/pull/470>`_.

.. rubric:: Bugfixes

- Fix a bug when computing the Schur decomposition for reducible Markov chains
  (*Schur vectors appear to not be D-orthogonal*). GPCCA requires the leading Schur vectors to be orthogonal w.r.t. a
  symmetric, positive definite matrix :math:`D` `PR 453 <https://github.com/theislab/cellrank/pull/453>`_.
- Fix not falling back to ``mode='monte_carlo'`` if no :mod:`jax` is found when using ``mode='stochastic'`` in
  :meth:`cellrank.tl.kernels.VelocityKernel.compute_transition_matrix`
  `PR 472 <https://github.com/theislab/cellrank/pull/472>`_.
- Fix :mod:`pandas` ``v1.0.1`` indexing error in :func:`cellrank.tl.lineage_drivers`
  `PR 475 <https://github.com/theislab/cellrank/pull/475>`_.
- Fix not correctly propagating colors during aggregation in :class:`cellrank.tl.Lineage`
  `PR 482 <https://github.com/theislab/cellrank/pull/482>`_.

1.1.0 :small:`2020-11-17`
~~~~~~~~~~~~~~~~~~~~~~~~~
This release includes:

.. rubric:: Additions

- :func:`cellrank.tl.lineage_drivers` computes p-values for the identified driver genes now, using either
  a Fisher-transformation to approximate the distribution of the test statistic under the null hypothesis
  or an exact, permutation based test. Corrects for multiple-testing.
- :meth:`cellrank.tl.kernels.VelocityKernel.compute_transition_matrix` now allows different metrics to be used to
  compare velocity vectors with expression-differences across neighboring cells. We add cosine-correlation and
  dot-product schemes and we allow the user to input their own scheme. It has been shown recently by :cite:`li:20`
  that the choice of metric can lead to slightly different results. Users can now also supply their own scheme as long
  as it follows the signature of :class:`cellrank.tl.kernels.SimilaritySchemeABC`.
- :func:`cellrank.datasets.reprogramming` has been added to allow for easy reproducibility of the time & memory
  benchmarking results in our `CellRank preprint <https://doi.org/10.1101/2020.10.19.345983>`_. This is a reprogramming
  dataset from :cite:`morris:18`.

.. rubric:: Bugfixes

- Fix not vendorizing correct :mod:`msmtools` which sometimes caused densification of a sparse matrix.
- Bump scanpy version requirement to 1.6 to fix plotting `PR 444 <https://github.com/theislab/cellrank/pull/444>`_.


1.0.0 :small:`2020-10-17`
~~~~~~~~~~~~~~~~~~~~~~~~~
- Fix a bug when subsetting :class:`cellrank.tl.Lineage`
- Add renaming terminal states :meth:`cellrank.tl.estimators.BaseEstimator.rename_terminal_states`
- Enable negative binomial distribution for :class:`cellrank.ul.models.GAMR`
- Remove previously deprecated functions
- Add :class:`cellrank.ul.models.FailedModel` inspired by the maybe monad
- Allow returning models when doing bulk fitting
- Add ``transpose`` parameter for :func:`cellrank.pl.gene_trends`
- Various minor bugfixes

1.0.0-rc.11 :small:`2020-09-25`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Rename ``final states`` to ``terminal states``
- Fix pickling if :class:`cellrank.tl.estimators.BaseEstimator`
- Fix various color bugs
- Improve :class:`cellrank.tl.kernels.PrecomputedKernel`
- Update gallery
- Other various minor changes

1.0.0-rc.0 :small:`2020-07-15`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Fix pickling of :class:`cellrank.tl.Lineage`
- Add additional options to :func:`cellrank.pl.heatmap`
- Updated documentation

1.0.0-b.8 :small:`2020-07-12`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Add installation options for PETSc and SLEPc
- Add iterative solver for absorption probabilities
- Add minor :class:`cellrank.tl.Lineage` improvements
- Fix docstring issues

1.0.0-b.2 :small:`2020-07-02`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Fix installation by including future-fstrings

1.0.0-b.1 :small:`2020-07-02`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Initial beta pre-release
