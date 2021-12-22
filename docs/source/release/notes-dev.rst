CellRank dev (2021-12-02)
=========================

Features
--------

- Update CytoTRACE calculation to check inf vals
  `#761 <https://github.com/theislab/cellrank/pull/761>`__

- Refactor :class:`cellrank.tl.kernels.TransportMapKernel` for easier definition of extensions.
  `#764 <https://github.com/theislab/cellrank/pull/764>`__

- Allow plotting negatively correlated genes in :func:`cellrank.pl.lineage_drivers`.
  `#779 <https://github.com/theislab/cellrank/pull/779>`__


Bugfixes
--------

- Fix :class:`cellrank.external.kernels.WOTKernel` saving OTModel and being unable to pickle itself.
  `#748 <https://github.com/theislab/cellrank/pull/748>`__

- Fix when associating term/macro states assignment when reference is not categorical string.
  `#750 <https://github.com/theislab/cellrank/pull/750>`__

- Fix color creating in :class:`cellrank.tl.kernels.ExperimentalTimeKernel`.
  `#784 <https://github.com/theislab/cellrank/pull/784>`__
