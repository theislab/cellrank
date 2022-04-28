CellRank dev (2022-04-28)
=========================

Features
--------

- Update CytoTRACE calculation to check inf vals
  `#761 <https://github.com/theislab/cellrank/pull/761>`__

- Refactor :class:`cellrank.tl.kernels.TransportMapKernel` for easier definition of extensions.
  `#764 <https://github.com/theislab/cellrank/pull/764>`__

- Allow plotting negatively correlated genes in :func:`cellrank.pl.lineage_drivers`.
  `#779 <https://github.com/theislab/cellrank/pull/779>`__

- Refactor kernel module. These changes include:
  - allow kernel elementwise multiplication
  - decouple CytoTRACE computation from kernel initialization, see :meth:`cellrank.tl.kernels.CytoTRACEKernel.compute_cytottrace`
  - change argument ``mode`` to ``model`` in :class:`cellrank.tl.kernels.VelocityKernel`
  - change argument ``scheme`` to ``similarity`` in :class:`cellrank.tl.kernels.VelocityKernel`
  - remove ``'sampling'`` option in :class:`cellrank.tl.kernels.VelocityKernel`
  - distinguish between unidirectional and bidirectional kernels
  - bidirectional kernel inversion is no longer inplace, returns a new object
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


Miscellaneous
-------------

- Ignore :mod:`pygam` deprecation warnings.
  `#798 <https://github.com/theislab/cellrank/pull/798>`__

- Fix PETSc/SLEPc in CI for Linux. For now, macOS PETSc/SLEPc remains disabled.
  `#850 <https://github.com/theislab/cellrank/pull/850>`__
