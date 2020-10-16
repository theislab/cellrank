Release Notes
=============

.. role:: small

Version 1.0
-----------

1.0.0 :small:`2020-10-17`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Fix bug when subsetting :class:`cellrank.tl.Lineage`
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
