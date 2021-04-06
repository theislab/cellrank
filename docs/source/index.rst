|PyPI| |Bioconda| |Downloads| |CI| |Notebooks| |Docs| |Codecov|

CellRank for directed single-cell fate mapping
========================================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/cellrank_fate_map.png
   :width: 600px
   :align: center

**CellRank** is a toolkit to uncover cellular dynamics based on Markov state modelling of single-cell data,
taking into account the stochastic nature of cellular fate decisions. It contains two main modules:
**kernels** compute cell-cell transition probabilities and **estimators** generate hypothesis based on these
transition probabilities. Our kernels work with a variety of input data including RNA velocity (see [Manno18]_
and [Bergen20]_), transcriptomic similarity, spatial proximity and pseudotime. Our ``VelocityKernel`` takes into
account **uncertainty in the velocities** and allows you to aggregate the short-range fate relations given by RNA
velocity into longer trends: from initial to terminal states along the phenotypic manifold. Our main estimator
is *Generalized Perron Cluster Cluster Analysis* (G-PCCA) [GPCCA18]_, implemented in the novel `pyGPCCA`_ package.
GPCCA coarse-grains the Markov chain into a set of macrostates which represent initial and terminal states,
as well as transient intermediate states. For each transient cell, i.e. for each cell that's not assigned to a
terminal state, we then compute its fate probability of it reaching any of the terminal states.
We show an example of such a fate map in the figure above, which has been computed using the data
of [Panc19]_.

CellRank scales to large cell numbers, is fully compatible with `scanpy`_ and `scvelo <https://scvelo.readthedocs.io/>`_
and is easy to use. To get started, see our `tutorial`_.

Manuscript
^^^^^^^^^^
Please see our `preprint`_ on **bioRxiv** to learn more.

Latest additions
^^^^^^^^^^^^^^^^

.. include:: latest_additions.rst

CellRank's key applications
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- compute initial & terminal as well as intermediate macrostates of your biological system
- infer fate probabilities towards the terminal states for each individual cell
- visualize gene expression trends along specific lineages while accounting for the continuous nature of
  fate determination
- identify potential driver genes for each identified cellular trajectory

Why is it called "CellRank"?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CellRank **does not** rank cells, we gave the package this name because just like Google's original `PageRank`_
algorithm, it works with Markov chains to aggregate relationships between individual objects (cells vs. websites)
to learn about more global properties of the underlying dynamics (initial & terminal states and fate probabilities vs.
website relevance).

Support
^^^^^^^
We welcome your feedback! Feel free to open an `issue <https://github.com/theislab/cellrank/issues/new/choose>`_
or send us an `email <mailto:info@cellrank.org>`_ if you encounter a bug, need our help or just want to make a
comment/suggestion.

CellRank was developed in collaboration between the `Theislab`_ and the `Peerlab`_.

.. toctree::
    :caption: General
    :maxdepth: 2
    :hidden:

    installation
    api
    external_api
    classes
    release_notes
    references

.. toctree::
   :caption: Gallery
   :maxdepth: 2
   :hidden:

   auto_examples/index.rst

.. toctree::
   :caption: Tutorials
   :maxdepth: 2
   :hidden:

   cellrank_basics
   kernels_and_estimators
   beyond_rna_velocity
   creating_new_kernel

.. |PyPI| image:: https://img.shields.io/pypi/v/cellrank.svg
    :target: https://pypi.org/project/cellrank
    :alt: PyPI

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/cellrank
    :target: https://bioconda.github.io/recipes/cellrank/README.html
    :alt: Bioconda

.. |CI| image:: https://img.shields.io/github/workflow/status/theislab/cellrank/CI/master
    :target: https://github.com/theislab/cellrank/actions
    :alt: CI

.. |Notebooks| image:: https://img.shields.io/github/workflow/status/theislab/cellrank_notebooks/CI/master?label=notebooks
    :target: https://github.com/theislab/cellrank_notebooks/actions
    :alt: CI-Notebooks

.. |Docs|  image:: https://img.shields.io/readthedocs/cellrank
    :target: https://cellrank.readthedocs.io/en/stable
    :alt: Documentation

.. |Downloads| image:: https://pepy.tech/badge/cellrank
    :target: https://pepy.tech/project/cellrank
    :alt: Downloads

.. |Codecov| image:: https://codecov.io/gh/theislab/cellrank/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/theislab/cellrank
    :alt: Coverage

.. _preprint: https://www.biorxiv.org/content/10.1101/2020.10.19.345983v1
.. _PageRank: https://en.wikipedia.org/wiki/PageRank#cite_note-1
.. _tutorial: https://cellrank.readthedocs.io/en/stable/cellrank_basics.html
.. _scanpy: https://scanpy.readthedocs.io/en/stable/
.. _Theislab: https://www.helmholtz-muenchen.de/icb/research/groups/theis-lab/overview/index.html
.. _Peerlab: https://www.mskcc.org/research/ski/labs/dana-pe-er
.. _pyGPCCA: https://pygpcca.readthedocs.io/en/latest/
