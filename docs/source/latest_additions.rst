.. role:: small

1.1.0 :small:`2020-11-14`
~~~~~~~~~~~~~~~~~~~~~~~~~

This release includes TODO:

.. rubric:: Bugfixes

- Fix msmtools version
- scanpy 1.6 requirement to fix a plotting FN

.. rubric:: Additions

- :func:`cellrank.tl.lineage_drivers` computes p-values for the identified driver genes now, using either
a Fisher-transformation to approximate the distribution of the test statistic under the null hypothesis
or an exact, permutation based test. Corrects for multiple-testing.
- :class:`cellrank.tl.kernels.VelocityKernel` now allows different metrics to be used to compare velocity
vectors with expression-differences across neighboring cells. We add cosine-correlation and dot-product
schemes and we allow the user to input their own scheme. It has been shown recently by `Li et al. <https://doi.org/10.1101/2020.09.19.304584>`_
that the choice of metric can lead to slightly different results.
- :func:`cellrank.datasets.reprogramming` has been added to allow for easy reproducibility of the time & memory
 benchmarking results in our `CellRank preprint <https://doi.org/10.1101/2020.10.19.345983>`_. This is a reprogramming
dataset from [Morris18]_.
