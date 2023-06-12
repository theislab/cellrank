CellRank 1.1.0 (2020-11-17)
===========================

Bugfixes
--------

- Fix not vendorizing correct :mod:`msmtools` which sometimes densified sparse matrix.
- Bump scanpy version requirement to 1.6 to fix plotting `PR 444 <https://github.com/theislab/cellrank/pull/444>`_.

Features
--------

- :func:`cellrank.tl.lineage_drivers` computes p-values for the identified driver genes now, using either
  a Fisher-transformation to approximate the distribution of the test statistic under the null hypothesis
  or an exact, permutation based test. Corrects for multiple-testing.
- :meth:`cellrank.tl.kernels.VelocityKernel.compute_transition_matrix` now allows different metrics to be used to
  compare velocity vectors with expression-differences across neighboring cells. We add cosine-correlation and
  dot-product schemes and we allow the user to input their own scheme. It has been shown recently by :cite:`li:20`
  that the choice of metric can lead to slightly different results. Users can now also supply their own scheme as long
  as it follows the signature of :class:`cellrank.tl.kernels.SimilaritySchemeABC`.
- :func:`cellrank.datasets.reprogramming` has been added to allow for easy reproducibility of the time & memory
  benchmarking results in our manuscript :cite:`lange:22`. This is a reprogramming dataset from :cite:`morris:18`.
