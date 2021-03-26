.. role:: small

1.3.0 :small:`2021-03-29`
~~~~~~~~~~~~~~~~~~~~~~~~~

This release includes some major additions which make CellRank more applicable with and without RNA velocity
information. In particular, it includes:

.. rubric:: Additions

- Add new kernel :class:`cellrank.tl.kernels.CytoTRACEKernel` which computes cell-cell transition probabilities based
  on the CytoTRACE score [Cyto20]_, a measure of differentiation potential,
  `PR 527 <https://github.com/theislab/cellrank/pull/527>`_.
- Add external API :mod:`cellrank.external` with a stationary optimal transport kernel
  :class:`cellrank.external.kernels.OTKernel` contributed from [Zhang21]_, as well as a
  `contributing guide <https://github.com/theislab/cellrank/blob/master/CONTRIBUTING.rst>`_,
  `PR 522 <https://github.com/theislab/cellrank/pull/522>`_.
- Rename ``cellrank.tl.kernels.PalantirKernel`` to :class:`cellrank.tl.kernels.PseudotimeKernel` and add
  hard thresholding scheme inspired by [Setty19]_, a soft thresholding scheme inspired by [VIA21]_ and a custom scheme
  when computing the transition matrix, see e.g. :class:`cellrank.tl.kernels.SoftThresholdScheme`
  `PR 514 <https://github.com/theislab/cellrank/pull/514>`_.
- Add more flexibility to :class:`cellrank.tl.kernels.ConnectivityKernel`, allowing it to use any cell-cell similarities
  from :attr:`anndata.AnnData.obsp`, such as spatial similarities from :mod:`squidpy` [Palla21]_
  `PR 501 <https://github.com/theislab/cellrank/pull/501>`_.
- Revamp `Pancreas Advanced <https://cellrank.readthedocs.io/en/latest/pancreas_advanced.html>`_ tutorial
  to showcase CellRank's modular structure of kernels and estimators.
  `PR 32 <https://github.com/theislab/cellrank_notebooks/pull/32>`_.
- Add 2 new tutorials:

  - `Beyond RNA velocity <https://cellrank.readthedocs.io/en/latest/beyond_rna_velocity.html>`_: shows how to use
    CellRank when no RNA velocity information is available.
    `PR 32 <https://github.com/theislab/cellrank_notebooks/pull/32>`_
  - `Creating a new kernel <https://cellrank.readthedocs.io/en/latest/creating_new_kernel.html>`_: explains how to
    create your own custom kernel class that estimates cell-cell transition probabilities
    `PR 31 <https://github.com/theislab/cellrank_notebooks/pull/31>`_.

- Add projection of transition matrix onto an embedding :meth:`cellrank.tl.kernels.Kernel.compute_projection`
- Add random walk simulation and visualization in an embedding :meth:`cellrank.tl.kernels.Kernel.plot_random_walks`
  `PR 537 <https://github.com/theislab/cellrank/pull/537>`_.
- Add :meth:`cellrank.tl.Lineage.priming_degree` `PR 502 <https://github.com/theislab/cellrank/pull/502>`_
  which estimates a cell's plasticity/differentiation potential based on ideas by [Setty19]_ and [Velten17]_.
- Add checks for transition matrix irreducibility `PR 516 <https://github.com/theislab/cellrank/pull/516>`_.
- Add Zebrafish development dataset from [Farrel18]_ `PR 539 <https://github.com/theislab/cellrank/pull/539>`_.
- Speed-up stationary distribution calculation in :mod:`pygpcca` `PR 22 <https://github.com/msmdev/pyGPCCA/pull/22>`_.

.. rubric:: Bugfixes

- Fix various bugs when plotting multiple gene trends `PR 487 <https://github.com/theislab/cellrank/pull/487>`_.
- Fix gene trend smoothing not working for 1 lineage `PR 512 <https://github.com/theislab/cellrank/pull/512>`_.
- Fix :mod:`pandas` error when computing macrostates `PR 513 <https://github.com/theislab/cellrank/pull/513>`_.
- Remove malfunctioning *Edit on GitHub* from the documentation
  `PR 538 <https://github.com/theislab/cellrank/pull/538>`_.
