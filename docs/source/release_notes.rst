Release Notes
=============

.. role:: small

Version 1.0
-----------

1.1.0 :small:`2020-11-17`
~~~~~~~~~~~~~~~~~~~~~~~~~

This release includes:

.. rubric:: Bugfixes

- Fix not vendorizing correct :mod:`msmtools` which sometimes caused densification of a sparse matrix.
- Bump scanpy version requirement to 1.6 to fix plotting `PR 444 <https://github.com/theislab/cellrank/pull/444>`_.

.. rubric:: Additions

- :func:`cellrank.tl.lineage_drivers` computes p-values for the identified driver genes now, using either
  a Fisher-transformation to approximate the distribution of the test statistic under the null hypothesis
  or an exact, permutation based test. Corrects for multiple-testing.
- :meth:`cellrank.tl.kernels.VelocityKernel.compute_transition_matrix` now allows different metrics to be used to
  compare velocity vectors with expression-differences across neighboring cells. We add cosine-correlation and
  dot-product schemes and we allow the user to input their own scheme. It has been shown recently by [Li2020]_
  that the choice of metric can lead to slightly different results. Users can now also supply their own scheme as long
  as it follows the signature of :class:`cellrank.tl.kernels.SimilaritySchemeABC`.
- :func:`cellrank.datasets.reprogramming` has been added to allow for easy reproducibility of the time & memory
  benchmarking results in our `CellRank preprint <https://doi.org/10.1101/2020.10.19.345983>`_. This is a reprogramming
  dataset from [Morris18]_.

1.0.0 :small:`2020-10-17`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Fix a bug when subsetting :class:`cellrank.tl.Lineage`
- Add renaming terminal states :meth:`cellrank.tl.estimators.BaseEstimator.rename_terminal_states`
- Enable negative binomial distribution for :class:`cellrank.ul.models.GAMR`
- Remove previously deprecated functions
- Add :class:`cellrank.ul.models.FailedModel` inspired by the maybe monad
- Allow returning models when doing bulk fitting
- Add ``transpose`` parameter for :func:`cellrank.pl.gene_trends`
- Various other minor bugfixes

1.0.0-rc.11 :small:`2020-09-25`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Rename ``metastable states`` to ``macrostates``
- Rename ``final states`` to ``terminal states``
- Fix pickling if :class:`cellrank.tl.estimators.BaseEstimator`
- Fix various color bugs
- Improve :class:`cellrank.tl.kernels.PrecomputedKernel`
- Update gallery
- Other various minor changes

1.0.0-rc.0 :small:`2020-07-15`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Fix pickling of :class:`cellrank.tl.Lineage` improvements
- Add additional options to :func:`cellrank.pl.heatmap`
- Updated documentation

1.0.0-b.8 :small:`2020-07-12`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Add installation options for PETSc and SLEPc
- Add iterative solver for absorption proabilities
- Add minor :class:`cellrank.tl.Lineage` improvements
- Fix docstring issues

1.0.0-b.2 :small:`2020-07-02`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Fix installation by including future-fstrings.

1.0.0-b.1 :small:`2020-07-02`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Initial beta pre-release.
