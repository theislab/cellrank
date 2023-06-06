CellRank 1.0.0 (2020-10-17)
===========================

Bugfixes
--------

- Fix a bug when subsetting :class:`cellrank.tl.Lineage`
- Remove previously deprecated functions

Features
--------

- Add renaming terminal states :meth:`cellrank.tl.estimators.BaseEstimator.rename_terminal_states`
- Enable negative binomial distribution for :class:`cellrank.ul.models.GAMR`
- Add :class:`cellrank.ul.models.FailedModel` inspired by the maybe monad
- Allow returning models when doing bulk fitting
- Add ``transpose`` parameter for :func:`cellrank.pl.gene_trends`
