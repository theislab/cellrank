.. role:: small

1.2.0 :small:`2021-02-02`
~~~~~~~~~~~~~~~~~~~~~~~~~

This release includes:

.. rubric:: Bugfixes

- Fix a bug when computing the Schur decomposition for reducible Markov chains
  (*Schur vectors appear to not be D-orthogonal*). GPCCA requires the leading Schur vectors to be orthogonal w.r.t. a
  symmetric, positive definite matrix :math:`D` `PR 453 <https://github.com/theislab/cellrank/pull/453>`_.
- Fix not falling back to ``mode='monte_carlo'`` if no :mod:`jax` is found when using ``mode='stochastic'`` in
  :meth:`cellrank.tl.kernels.VelocityKernel.compute_transition_matrix`
  `PR 472 <https://github.com/theislab/cellrank/pull/472>`_.
- Fix :mod:`pandas` ``v1.0.1`` indexing error in :func:`cellrank.tl.lineage_drivers`
  `PR 475 <https://github.com/theislab/cellrank/pull/475>`_.
- Fix not correctly propagating colors during aggregation in :class:`cellrank.tl.Lineage`
  `PR 482 <https://github.com/theislab/cellrank/pull/482>`_.

.. rubric:: Additions

- Completely **refactored the underlying code base of GPCCA** and set it up as it's own package called
  `pyGPCCA <https://pygpcca.readthedocs.io/en/latest/>`_ with documentation and an example. Going forwards, this will
  ensure that one of the "engines" of CellRank is also easy to maintain to extend. Further, this will make CellRank's
  installation more convenient by not needing to vendorize additional dependencies
  `PR 472 <https://github.com/theislab/cellrank/pull/472>`_.
- Add :func:`cellrank.pl.circular_projection` visualizing computed fate probabilities as done in [Velten17]_,
  see :ref:`sphx_glr_auto_examples_plotting_plot_circular_embedding.py`.
  `PR 459 <https://github.com/theislab/cellrank/pull/459>`_.
- Allow legends not to be plotted by passing ``legend_loc="none"``, as done in `scVelo <https://scvelo.org>`_
  `PR 470 <https://github.com/theislab/cellrank/pull/470>`_.
