CellRank dev (2021-09-13)
=========================

Features
--------

- Add ``threshold`` to :meth:`cellrank.external.kernels.WOTKernel.compute_transition_matrix` to increase sparsity and
  speed up :class:`cellrank.estimators.GPCCA`.
  `#696 <https://github.com/theislab/cellrank/pull/696>`__

- Allow using a column from :attr:`anndata.AnnData.var` as gene symbols for some plotting functions.
  `#726 <https://github.com/theislab/cellrank/pull/726>`__


Bugfixes
--------

- Fix towncrier release generation CI.
  `#701 <https://github.com/theislab/cellrank/pull/701>`__

- Update ``towncrier`` to display development release notes.
  `#709 <https://github.com/theislab/cellrank/pull/709>`__

- Restricts computation of embedding projection to kNN based kernels.
  `#733 <https://github.com/theislab/cellrank/pull/733>`__

- Fix :meth:`cellrank.external.kernels.WOTKernel.compute_transition_matrix` silently ignoring unexpected kwargs.
  `#737 <https://github.com/theislab/cellrank/pull/737>`__

- Use actual number of nearest neighbors in :class:`cellrank.kernels.PseudotimeKernel`
  when using hard thresholding scheme.
  `#738 <https://github.com/theislab/cellrank/pull/738>`__


Deprecations (in next major release)
------------------------------------

- Deprecate :mod:`cellrank.tl`, including the high level API and rename
  :mod:`cellrank.ul.models` to :mod:`cellrank.models`.
  `#695 <https://github.com/theislab/cellrank/pull/695>`__


Miscellaneous
-------------

- Fix many test warnings.
  `#704 <https://github.com/theislab/cellrank/pull/704>`__

- Speed-up testing by not using stochastic mode in :class:`cellrank.kernels.VelocityKernel` where not necessary.
  `#705 <https://github.com/theislab/cellrank/pull/705>`__

- Enable ``tox`` in CI.
  `#713 <https://github.com/theislab/cellrank/pull/713>`__

- Update deployment CI and CONTRIBUTING.rst based on a new branching structure.
  `#725 <https://github.com/theislab/cellrank/pull/725>`__

- Add Python 3.9 CI testing.
  `#730 <https://github.com/theislab/cellrank/pull/730>`__


Documentation
-------------

- Add Google Colab links for tutorials.
  `#707 <https://github.com/theislab/cellrank/pull/707>`__

- Allow ``towncrier`` to generate bleeding-edge development notes.
  `#712 <https://github.com/theislab/cellrank/pull/712>`__
