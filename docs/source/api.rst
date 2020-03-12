API
===
CellRank offers two possibilities to interact with its core functionality: either
directly through the :class:`cellrank.tl.MarkovChain`
and :class:`cellrank.tl.kernels.Kernel` classes, or through the high
level functions :func:`cellrank.tl.root_final` and :func:`cellrank.tl.lineages`.
The former option offers full control whereas the latter is simpler to handle for new users.

Additionally, there is a set of plotting functions which can be used, no matter
which interface was chosen to uncover the start-/endpoints and lineages. See our
tutorials to learn more.

The utilities are mainly for fitting continuous models to gene expression data
and are utilised in some of the plotting functions, like :func:`cellrank.pl.gene_trends`.

Tools
~~~~~

.. module:: cellrank.tl
.. currentmodule:: cellrank

.. autosummary::

    :toctree: api/tl

    tl.MarkovChain
    tl.partition
    tl.root_final
    tl.lineages
    tl.gene_importance
    tl.transition_matrix
    tl.kernels.VelocityKernel
    tl.kernels.ConnectivityKernel
    tl.kernels.PalantirKernel

Plotting
~~~~~~~~

.. module:: cellrank.pl
.. currentmodule:: cellrank

.. autosummary::

    :toctree: api/pl

    pl.gene_trends
    pl.heatmap
    pl.cluster_lineage
    pl.cluster_fates
    pl.similarity_plot
    pl.graph
    pl.lineages
    pl.composition

Utilities
~~~~~~~~~

.. module:: cellrank.ul
.. currentmodule:: cellrank

.. autosummary::

    :toctree: api/ul

    ul.models.SKLearnModel
    ul.models.GamMGCVModel
