CellRank dev (2023-06-08)
=========================

Features
--------

- Update CytoTRACE calculation to check infinite values.
  `#761 <https://github.com/theislab/cellrank/pull/761>`__

- Refactor :class:`cellrank.tl.kernels.TransportMapKernel` for easier definition of extensions.
  `#764 <https://github.com/theislab/cellrank/pull/764>`__

- Allow plotting negatively correlated genes in :func:`cellrank.pl.lineage_drivers`.
  `#779 <https://github.com/theislab/cellrank/pull/779>`__

- Refactor kernel module. These changes include:
  - allow kernel element-wise multiplication
  - decouple CytoTRACE computation from kernel initialization, see :meth:`cellrank.tl.kernels.CytoTRACEKernel.compute_cytottrace`
  - change argument ``mode`` to ``model`` in :class:`cellrank.tl.kernels.VelocityKernel`
  - change argument ``scheme`` to ``similarity`` in :class:`cellrank.tl.kernels.VelocityKernel`
  - remove ``'sampling'`` option in :class:`cellrank.tl.kernels.VelocityKernel`
  - distinguish between unidirectional and bidirectional kernels
  - bidirectional kernel inversion is no longer in-place, returns a new object
  - add more memory efficient way of copying kernels
  - improve :class:`cellrank.kernels.PrecomputedKernel`
  - abstract kernel-related plotting functions (projection and random walks)
  - temporarily remove kernel tricks examples
  - fix various typos in docstrings
  - fix corner-case bug when computing absorption probabilities
  - fix not setting transport maps in :class:`cellrank.tl.kernels.TransportMapKernel`
  - fix automatic threshold removing 1 row in :class:`cellrank.tl.kernels.TransportMapKernel`
  `#791 <https://github.com/theislab/cellrank/pull/791>`__

- Allow local block-diagonal connectivities in :class:`cellrank.tl.kernels.TransportMapKernel`. Rename ``last_time_point`` argument to ``self_transitions`` and add ``conn_weight``. ``self_transitions`` can also be applied to specific diagonal blocks only.
  `#828 <https://github.com/theislab/cellrank/pull/828>`__

- Add :func:`cellrank.datasets.bone_marrow`.
  `#886 <https://github.com/theislab/cellrank/pull/886>`__

- Add argument ``subset_to_serum`` to :func:`cellrank.datasets.reprogramming_schiebinger` to allow downloading the subsetted data. This includes the transition matrix computed with the :class:`cellrank.external.WOTKernel`.
  `#890 <https://github.com/theislab/cellrank/pull/890>`__

- Updates kernels/estimators ``__str__/__repr__`` and allows for better method chaining.
  `#896 <https://github.com/theislab/cellrank/pull/896>`__

- Allow passing connectivities for transition matrix projection. Useful when the kernel is not kNN-based.
  `#930 <https://github.com/theislab/cellrank/pull/930>`__


Bugfixes
--------

- Fix :class:`cellrank.external.kernels.WOTKernel` saving OTModel and being unable to pickle itself.
  `#748 <https://github.com/theislab/cellrank/pull/748>`__

- Fix when associating term/macro states assignment when reference is not categorical string.
  `#750 <https://github.com/theislab/cellrank/pull/750>`__

- Fix color creating in :class:`cellrank.tl.kernels.ExperimentalTimeKernel`.
  `#784 <https://github.com/theislab/cellrank/pull/784>`__

- Fix all q-values being NaN if 1 p-value was NaN. Also warn and set to NaN instead of raise when correlations are not in ``[0, 1]]`` interval.
  `#835 <https://github.com/theislab/cellrank/pull/835>`__

- Do not save enums in :class:`anndata.AnnData`.
  `#842 <https://github.com/theislab/cellrank/pull/842>`__

- Fix :class:`cellrank.tl.Lineage` subsetting (non-existent overlapping keys) and refactor the implementation.
  `#861 <https://github.com/theislab/cellrank/pull/861>`__

- Adds references to the docs
  `#887 <https://github.com/theislab/cellrank/pull/887>`__

- Fix computing time to absorption to aggregated terminal states.
  `#923 <https://github.com/theislab/cellrank/pull/923>`__

- Fix initial/terminal points size in random walk.
  `#929 <https://github.com/theislab/cellrank/pull/929>`__

- This adapts some of our docs for the 2.0 release; it clears up the landing page and introduces some new pages, including the "teams", "how to cite us", and "about cellrank" pages.
  `#936 <https://github.com/theislab/cellrank/pull/936>`__

- This just adds a new concept figure to be used in the updated docs.
  `#937 <https://github.com/theislab/cellrank/pull/937>`__

- Update to the overview figure
  `#938 <https://github.com/theislab/cellrank/pull/938>`__

- Import ``Iterable`` from ``collections.abc`` and not ``collections``. ``Iterable`` was removed from ``collections`` in
  Python 3.10
  `#943 <https://github.com/theislab/cellrank/pull/943>`__

- Require ``Python >= 3.8`` and fix :func:`cellrank.pl.log_odds` colormap type.
  `#949 <https://github.com/theislab/cellrank/pull/949>`__

- Allow NaN values when determining foreground/background color.
  `#951 <https://github.com/theislab/cellrank/pull/951>`__

- Fix sparse array conversion in ``networkx>=3.0``.
  `#978 <https://github.com/theislab/cellrank/pull/978>`__

- Removes default value for ``time_key``.
  `#989 <https://github.com/theislab/cellrank/pull/989>`__

- Fix cluster subsetting in :func:`cellrank.pl.aggregate_absorption_probabilities` and ``mode = 'violin'``.
  `#1007 <https://github.com/theislab/cellrank/pull/1007>`__

- Fix not being able to save :class:`~anndata.AnnData` after :class:`~cellrank.kernels.CytoTRACEKernel` was run.
  `#1019 <https://github.com/theislab/cellrank/pull/1019>`__

- Fix not passing ``allow_overlap`` in some places.
  `#1020 <https://github.com/theislab/cellrank/pull/1020>`__

- Very minor fix to adapt the allowed options in macrostate plotting.
  `#1029 <https://github.com/theislab/cellrank/pull/1029>`__

- Fix ``adjustText`` modifying data in-place.
  `#1034 <https://github.com/theislab/cellrank/pull/1034>`__

- Adds a reference for Hhex
  `#1035 <https://github.com/theislab/cellrank/pull/1035>`__

- Fix lineage color map in :func:`cellrank.pl.gene_trends`.
  `#1043 <https://github.com/theislab/cellrank/pull/1043>`__


Miscellaneous
-------------

- Ignore :mod:`pygam` deprecation warnings.
  `#798 <https://github.com/theislab/cellrank/pull/798>`__

- Fix PETSc/SLEPc in CI for Linux. For now, macOS PETSc/SLEPc remains disabled.
  `#850 <https://github.com/theislab/cellrank/pull/850>`__

- Remove default real-time and pseudotime arguments from :class:`cellrank.kernels.ExperimentalTimeKernel` and :class:`cellrank.kernels.PseudotimeKernel`.
  `#894 <https://github.com/theislab/cellrank/pull/894>`__

- Fix not being able to use ``minChi`` in :meth:`cellrank.estimators.GPCCA.fit`. Also compute 20 Schur vectors by default.
  `#913 <https://github.com/theislab/cellrank/pull/913>`__

- Prefer plotting macrostates/terminal in a discrete way.
  `#914 <https://github.com/theislab/cellrank/pull/914>`__

- Use ``", "`` instead of ``" or "`` when joining states.
  `#954 <https://github.com/theislab/cellrank/pull/954>`__

- Use :mod:`rpy2`'s local converter.
  `#1008 <https://github.com/theislab/cellrank/pull/1008>`__

- TODO.
  `#1040 <https://github.com/theislab/cellrank/pull/1040>`__


Documentation
-------------

- Add *Discourse* badge for questions and discussions.
  `#880 <https://github.com/theislab/cellrank/pull/880>`__

- Adds references to the ``bibtex`` file.
  `#883 <https://github.com/theislab/cellrank/pull/883>`__

- Add ``furo`` theme based on ``scvi-tools``.
  `#885 <https://github.com/theislab/cellrank/pull/885>`__

- Fix ``RTD`` build - submodules.
  `#901 <https://github.com/theislab/cellrank/pull/901>`__

- Fix CI badges and tox.
  `#975 <https://github.com/theislab/cellrank/pull/975>`__

- Polish team page.
  `#976 <https://github.com/theislab/cellrank/pull/976>`__
