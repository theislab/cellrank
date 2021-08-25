CellRank dev (2021-08-21)
=========================

Features
--------

- Allow using a column from :attr:`anndata.AnnData.var` as gene symbols for some plotting functions.
  `#726 <https://github.com/theislab/cellrank/pull/726>`__


Bugfixes
--------

- Fix towncrier release generation CI.
  `#701 <https://github.com/theislab/cellrank/pull/701>`__

- Update towncrier to display development release notes.
  `#709 <https://github.com/theislab/cellrank/pull/709>`__


Deprecations (in next major release)
------------------------------------

- Deprecate :mod:`cellrank.tl`, including the high level API and
  rename :mod:`cellrank.ul.models` to :mod:`cellrank.models`.
  `#695 <https://github.com/theislab/cellrank/pull/695>`__


Miscellaneous
-------------

- Fix many test warnings.
  `#704 <https://github.com/theislab/cellrank/pull/704>`__

- Speed-up testing by not using stochastic mode in :class:`cellrank.tl.kernels.VelocityKernel` where not necessary.
  `#705 <https://github.com/theislab/cellrank/pull/705>`__

- Enable ``tox`` in CI.
  `#713 <https://github.com/theislab/cellrank/pull/713>`__


Documentation
-------------

- Add Google Colab links for tutorials.
  `#707 <https://github.com/theislab/cellrank/pull/707>`__

- Allow towncrier to generate bleeding-edge development notes.
  `#712 <https://github.com/theislab/cellrank/pull/712>`__
