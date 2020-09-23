|PyPI| |Bioconda| |Downloads| |Travis| |Notebooks| |Docs| |Codecov|


CellRank - Probabilistic Fate Mapping using RNA Velocity
========================================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/cellrank_fate_map.png
   :width: 600px
   :align: center

**CellRank** is a toolkit to uncover cellular dynamics based on scRNA-seq data with RNA velocity annotation,
see `La Manno et al. (2018)`_ and `Bergen et al. (2020)`_. CellRank models cellular dynamics as a Markov chain, where transition
probabilities are computed based on RNA velocity and transcriptomic similarity, taking into account uncertainty
in the velocities. The Markov chain is coarse grained into a set of metastable states which represent root &
final states as well as transient intermediate states. For each cell, we obtain the probability of it belonging
to each metastable state, i.e. we compute a fate map on the single cell level. We show an example of such a fate
map in the figure above, which has been computed using the data of `pancreatic endocrinogenesis`_.

CellRank scales to large cell numbers, is fully compatible with `scanpy`_ and `scvelo`_ and is easy
to use. For **installation instructions**, **documentation** and **tutorials**, visit `cellrank.org`_.

CellRank's key applications
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- compute root & final as well as intermediate metastable states of your developmental/dynamical process
- infer fate probabilities towards these states for each single cell
- visualise gene expression trends towards/from specific states
- identify potential driver genes for each state

Installation
^^^^^^^^^^^^
Install CellRank by running::

    conda install -c conda-forge -c bioconda cellrank
    # or with extra libraries, useful for large datasets
    conda install -c conda-forge -c bioconda cellrank-krylov

or via PyPI::

    pip install cellrank
    # or with extra libraries, useful for large datasets
    pip install 'cellrank[krylov]'

Support
^^^^^^^
We welcome your feedback! Feel free to open an `issue <https://github.com/theislab/cellrank/issues/new/choose>`_
or send us an `email <mailto:info@cellrank.org>`_ if you encounter a bug, need our help or just want to make a
comment/suggestion.

CellRank was developed in collaboration between the `Theislab`_ and the `Peerlab`_.

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

.. _La Manno et al. (2018): https://doi.org/10.1038/s41586-018-0414-6

.. _Bergen et al. (2020): https://doi.org/10.1038/s41587-020-0591-3

.. _pancreatic endocrinogenesis: https://doi.org/10.1242/dev.173849

.. _scanpy: https://scanpy.readthedocs.io/en/latest/

.. _scvelo: https://scvelo.readthedocs.io/

.. _cellrank.org: http://cellrank.org

.. _Theislab: https://www.helmholtz-muenchen.de/icb/research/groups/theis-lab/overview/index.html

.. _Peerlab: https://www.mskcc.org/research/ski/labs/dana-pe-er
