|PyPI| |Bioconda| |Downloads| |CI| |Notebooks| |Docs| |Codecov|

CellRank - Probabilistic Fate Mapping using RNA Velocity
========================================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/cellrank_fate_map.png
   :width: 600px
   :align: center

**CellRank** is a toolkit to uncover cellular dynamics based on scRNA-seq data with RNA velocity annotation,
see `La Manno et al. (2018)`_ and `Bergen et al. (2020)`_. In short, CellRank models cellular dynamics as a
Markov chain, where transition probabilities are computed based on **RNA velocity and transcriptomic similarity**,
taking into account **uncertainty in the velocities** and the stochastic nature of cell fate decisions.
The Markov chain is coarse-grained into a set of macrostates which represent initial and terminal states,
as well as transient intermediate states using Generalized Perron Cluster Cluster Analysis (G-PCCA) [GPCCA18]_,
implemented in the novel `pyGPCCA`_ package. For each transient cell, i.e. for each cell that's not assigned to a
terminal state, we then compute its fate probability of it reaching any of the terminal states.
We show an example of such a fate map in the figure above, which has been computed using the data
of `pancreatic endocrinogenesis`_.

CellRank scales to **large cell numbers**, is fully compatible with `scanpy`_ and `scvelo`_ and is **easy to use**.
For **installation instructions**, **documentation** and **tutorials**, visit `cellrank.org`_.

Manuscript
^^^^^^^^^^
Please see our `preprint`_ on **bioRxiv** to learn more.

CellRank's key applications
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- compute **initial & terminal** as well as **intermediate** macrostates of your biological system
- infer **fate probabilities** towards the terminal states for each individual cell
- visualize **gene expression trends** along specific linegeages while accounting for the continous nature of fate
  determination
- identify **potential driver genes** for each identified cellular trajectory

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

Why is it called "CellRank"?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CellRank **does not** rank cells, we gave the package this name because just like Google's original `PageRank`_
algorithm, it works with Markov chains to aggregate relationships between individual objects (cells vs. websites)
to learn about more global properties of the underlying dynamics (initial & terminal states and fate probabilities vs.
website relevance).

Support
^^^^^^^
We welcome your feedback! Feel free to open an `issue <https://github.com/theislab/cellrank/issues/new/choose>`__
or send us an `email <mailto:info@cellrank.org>`_ if you encounter a bug, need our help or just want to make a
comment/suggestion.

CellRank was developed in collaboration between the `Theislab`_ and the `Peerlab`_.

.. |PyPI| image:: https://img.shields.io/pypi/v/cellrank.svg
    :target: https://pypi.org/project/cellrank
    :alt: PyPI

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/cellrank
    :target: https://bioconda.github.io/recipes/cellrank/README.html
    :alt: Bioconda

.. |CI| image:: https://img.shields.io/github/workflow/status/theislab/cellrank/CI/master
    :target: https://github.com/theislab/cellrank
    :alt: CI

.. |Notebooks| image:: https://img.shields.io/github/workflow/status/theislab/cellrank_notebooks/CI/master?label=notebooks
    :target: https://github.com/theislab/cellrank_notebooks
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

.. _preprint: https://doi.org/10.1101/2020.10.19.345983

.. _PageRank: https://en.wikipedia.org/wiki/PageRank#cite_note-1

.. _La Manno et al. (2018): https://doi.org/10.1038/s41586-018-0414-6

.. _Bergen et al. (2020): https://doi.org/10.1038/s41587-020-0591-3

.. _pancreatic endocrinogenesis: https://doi.org/10.1242/dev.173849

.. _scanpy: https://scanpy.readthedocs.io/en/latest/

.. _scvelo: https://scvelo.readthedocs.io/

.. _cellrank.org: https://cellrank.org

.. _Theislab: https://www.helmholtz-muenchen.de/icb/research/groups/theis-lab/overview/index.html

.. _Peerlab: https://www.mskcc.org/research/ski/labs/dana-pe-er

.. _pyGPCCA: https://pygpcca.readthedocs.io/en/latest/

.. _GPCCA18: https://doi.org/10.1021/acs.jctc.8b00079
