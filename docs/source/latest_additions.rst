.. role:: small

1.3.0 :small:`2021-03-23`
~~~~~~~~~~~~~~~~~~~~~~~~~

This release includes:

.. rubric:: Additions

- Add new kernel :class:`cellrank.tl.kernels.CytoTRACEKernel` which computes the transition probabilities based on a
  pseudotime derived from the [Cyto20]_ score `PR 527 <https://github.com/theislab/cellrank/pull/527>`_.
- Add external API :mod:`cellrank.external` with a stationary optimal transport kernel
  :class:`cellrank.external.kernels.OTKernel` contributed from [Zhang21]_, as well as a
  `contributing guide <https://github.com/theislab/cellrank/blob/master/CONTRIBUTING.rst>`_,
  `PR 522 <https://github.com/theislab/cellrank/pull/522>`_.
- Rename ``cellrank.tl.kernels.PalantirKernel`` to :class:`cellrank.tl.kernels.PseudotimeKernel` and add
  hard thresholding scheme mimicking [Setty19]_, a soft thresholding scheme mimicking [VIA21]_ and a custom scheme
  when computing the transition matrix, see e.g. :class:`cellrank.tl.kernels.SoftThresholdScheme`
  `PR 514 <https://github.com/theislab/cellrank/pull/514>`_.
- Add more flexibility to :class:`cellrank.tl.kernels.ConnectivityKernel`, allowing it to use any cell-cell similarities
  from :attr:`anndata.AnnData.obsp`, such as spatial similarities from :mod:`squidpy` [Palla21]_
  `PR 501 <https://github.com/theislab/cellrank/pull/501>`_.
- Revamp `Pancreas Advanced <https://cellrank.readthedocs.io/en/latest/pancreas_advanced.html>`_ tutorial
  to showcase TODO `PR TODO <TODO>`_.
- Add 2 new tutorials:

  - TODO
  - `Creating a new kernel <https://cellrank.readthedocs.io/en/latest/creating_new_kernel.html>`_: this short tutorial
    explains how to create your own custom kernel class that estimates cell-cell transition matrices
    `PR 31 <https://github.com/theislab/cellrank_notebooks/pull/31>`_.

- Add projection of transition matrix onto an embedding :meth:`cellrank.tl.kernels.Kernel.compute_projection`
- Add random walk simulation and visualization in an embedding :meth:`cellrank.tl.kernels.Kernel.plot_random_walks`
  `PR 537 <https://github.com/theislab/cellrank/pull/537>`_.
- Add :meth:`cellrank.tl.Lineage.priming_degree` `PR 502 <https://github.com/theislab/cellrank/pull/502>`_
  which calculates the cells' commitment based on the ideas from [Setty19]_ and [Velten17]_.
- Add checks for the transition matrix irreducibility `PR 516 <https://github.com/theislab/cellrank/pull/516>`_.

.. rubric:: Bugfixes

- Fix various bugs when plotting multiple gene trends `PR 487 <https://github.com/theislab/cellrank/pull/487>`_.
- Fix gene trend smoothing not working for 1 lineage `PR 512 <https://github.com/theislab/cellrank/pull/512>`_.
- Fix :mod:`pandas` error when computing macrostates `PR 513 <https://github.com/theislab/cellrank/pull/513>`_.
- Remove malfunctioning *Edit on GitHub* from the documentation
  `PR 538 <https://github.com/theislab/cellrank/pull/538>`_.
