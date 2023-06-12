CellRank 1.5.0 (2021-09-13)
===========================

Features
--------

- Add ``threshold`` to :meth:`cellrank.external.kernels.WOTKernel.compute_transition_matrix` to increase sparsity and
  speed up :class:`cellrank.tl.estimators.GPCCA`.
  `#696 <https://github.com/theislab/cellrank/pull/696>`__

- Allow using columns from :attr:`anndata.AnnData.var` as gene symbols for some plotting functions.
  `#726 <https://github.com/theislab/cellrank/pull/726>`__

- Completely refactor :mod:`cellrank.tl.estimators`. This includes the following changes:

      - define less coupled, more extensible class hierarchy.
      - follow more closely fit/predict paradigm of :mod:`sklearn`, see e.g.
        :meth:`cellrank.tl.estimators.GPCCA.fit` and :meth:`cellrank.tl.estimators.GPCCA.predict`.
      - make estimators implicitly maintain a more consistent state.
      - remove plotting of Schur vectors  and eigenvectors in an embedding.
      - remove cell-cycle warning for terminal states.
      - remove ``is_irreducible``, ``recurrent_classes`` and ``transient_classes`` properties.
      - remove optional irreducibility check from :meth:`cellrank.tl.estimators.GPCCA.compute_absorption_probabilities`
      - normalize estimator attribute names and key names when writing to :class:`anndata.AnnData`.
      - allow estimators to serialize self from/to :class:`anndata.AnnData`, see
        :meth:`cellrank.tl.estimators.BaseEstimator.from_adata` or :meth:`cellrank.tl.estimators.BaseEstimator.to_adata`.
      - allow estimators to be saved to a file without it :class:`anndata.AnnData`.
      - improve docstrings in various places.
      - write lineage drivers to :attr:`anndata.AnnData.varm` instead of :attr:`anndata.AnnData.var`.
      - fix various corner cases when solving linear systems (e.g. only 1 variable).

  In addition, add :meth:`cellrank.tl.Lineage.from_adata` to allow easy reconstruction of lineage objects.
  `#727 <https://github.com/theislab/cellrank/pull/727>`__


Bugfixes
--------

- Fix ``towncrier`` release generation CI.
  `#701 <https://github.com/theislab/cellrank/pull/701>`__

- Update ``towncrier`` to display development release notes.
  `#709 <https://github.com/theislab/cellrank/pull/709>`__

- Restricts computation of embedding projection to kNN based kernels.
  `#733 <https://github.com/theislab/cellrank/pull/733>`__

- Fix :meth:`cellrank.external.kernels.WOTKernel.compute_transition_matrix` silently ignoring unexpected kwargs.
  `#737 <https://github.com/theislab/cellrank/pull/737>`__

- Use actual number of nearest neighbors in :class:`cellrank.tl.kernels.PseudotimeKernel`
  when using hard threshold scheme.
  `#738 <https://github.com/theislab/cellrank/pull/738>`__

- Fix :func:`cellrank.pl.cluster_lineage` sometimes reusing the same ax.
  `#742 <https://github.com/theislab/cellrank/pull/742>`__


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

- Add *Google Colab* links for tutorials.
  `#707 <https://github.com/theislab/cellrank/pull/707>`__

- Allow ``towncrier`` to generate bleeding-edge development notes.
  `#712 <https://github.com/theislab/cellrank/pull/712>`__

- Fix docstrings and use typed returns.
  `#743 <https://github.com/theislab/cellrank/pull/743>`__
