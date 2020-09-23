|PyPI| |Bioconda| |Downloads| |Travis| |Notebooks| |Docs| |Codecov|


CellRank - Probabilistic Fate Mapping using RNA Velocity
===================================================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/cellrank_fate_map.png
   :width: 600px
   :align: center

**CellRank** is a toolkit to uncover cellular dynamics based on scRNA-seq data with RNA velocity annotation,
see [Manno18]_ and [Bergen20]_. CellRank models cellular dynamics as a Markov chain, where transition
probabilities are computed based on RNA velocity and transcriptomic similarity, taking into account uncertainty
in the velocities. The Markov chain is coarse grained into a set of metastable states which represent root &
final states as well as transient intermediate states. For each cell, we obtain the probability of it belonging
to each metastable state, i.e. we compute a fate map on the single cell level. We show an example of such a fate
map in the figure above, which has been computed using the data of [Panc19]_.

CellRank scales to large cell numbers, is fully compatible with `scanpy`_ and `scvelo`_ and is easy
to use. To get started, see our `tutorial`_.

CellRank's key applications
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- compute root & final as well as intermediate metastable states of your developmental/dynamical process
- infer fate probabilities towards these states for each single cell
- visualise gene expression trends towards/from specific states
- identify potential driver genes for each state

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

   pancreas_basic
   pancreas_advanced

.. |PyPI| image:: https://img.shields.io/pypi/v/cellrank.svg
    :target: https://pypi.org/project/cellrank
    :alt: PyPI

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/cellrank
    :target: https://bioconda.github.io/recipes/cellrank/README.html
    :alt: Bioconda

.. |Travis| image:: https://travis-ci.org/theislab/cellrank.svg?branch=master
    :target: https://travis-ci.com/github/theislab/cellrank
    :alt: CI

.. |Notebooks| image:: https://img.shields.io/travis/com/theislab/cellrank_notebooks?label=notebooks
    :target: https://travis-ci.com/github/theislab/cellrank_notebooks
    :alt: CI-Notebooks

.. |Docs|  image:: https://img.shields.io/readthedocs/cellrank
    :target: https://cellrank.readthedocs.io/en/latest
    :alt: Documentation

.. |Downloads| image:: https://pepy.tech/badge/cellrank
    :target: https://pepy.tech/project/cellrank
    :alt: Downloads

.. |Codecov| image:: https://codecov.io/gh/theislab/cellrank/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/theislab/cellrank
    :alt: Coverage

.. _tutorial: https://cellrank.readthedocs.io/en/latest/pancreas_basic.html

.. _PageRank: http://infolab.stanford.edu/~backrub/google.html

.. _documentation: https://cellrank.readthedocs.io

.. _scanpy: https://scanpy.readthedocs.io/en/latest/

.. _scvelo: https://scvelo.readthedocs.io/

.. _Theislab: https://www.helmholtz-muenchen.de/icb/research/groups/theis-lab/overview/index.html

.. _Peerlab: https://www.mskcc.org/research/ski/labs/dana-pe-er
