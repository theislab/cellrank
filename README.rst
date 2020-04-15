|Travis| |Docs| |Codecov| |License|


CellRank - Probabilistic Lineage Assignment using RNA Velocity
===================================================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/index_figure_endpoints.png
   :width: 600px
   :align: center

**CellRank** is a toolkit to uncover cellular development based on scRNA-seq data with RNA velocity annotation,
see `La Manno et al. (2018)`_ and `Bergen et al. (2019)`_.
In the figure above, we show the main features of CellRank applied to `pancreatic endocrinogenesis`_ -
starting from RNA velocities **(a)**, we infer root cells **(b)** and final cells **(c)**, and we compute
how likely each cell is to develop towards each of the identified groups of final cells **(d)**.
See our `tutorial`_ to learn how to apply CellRank to your own data.

CellRank utilizes the time derivative of gene expression given by RNA velocity to construct a Markov chain.
The information given by RNA velocity is combined with transcriptomic similarity and density corrected to yield
a robust estimate of cellular development directly in high dimensional gene expression space.
CellRank obtained its name due to conceptual similarities with `PageRank`_, Google’s original algorithm
for ranking web pages. Both algorithms construct a Markov Chain and use spectral methods to study its
long term evolution. Based on the Markov Chain, we infer root and final cells of development as well
as lineage probabilities, i.e. for each cell, the probability of it reaching any of the inferred
root and final cells of development. CellRank offers many possibilities to utilize these
lineage probabilities in downstream analysis, e.g. to investigate the fate of early cell clusters
or to plot gene expression trends along a given lineage.

CellRank requires Python3 version >= 3.6 and can be installed via::

    pip install git+https://github.com/theislab/cellrank

See the `documentation`_, which includes tutorial notebooks and a description of our complete API.

CellRank is fully compatible with `scanpy`_ and `scvelo`_ and was developed in collaboration
between the `Theislab`_ and the `Peerlab`_.


.. |Travis| image:: https://travis-ci.org/theislab/cellrank.svg?branch=master
    :target: https://travis-ci.org/theislab/cellrank

.. |Docs|  image:: https://img.shields.io/readthedocs/cellrank

.. |Codecov| image:: https://codecov.io/gh/theislab/cellrank/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/theislab/cellrank

.. |License| image:: https://img.shields.io/github/license/theislab/cellrank

.. _La Manno et al. (2018): https://doi.org/10.1038/s41586-018-0414-6

.. _Bergen et al. (2019): https://doi.org/10.1101/820936

.. _pancreatic endocrinogenesis: https://doi.org/10.1242/dev.173849

.. _tutorial: https://cellrank-notebooks.readthedocs.io/en/latest/pancreas_basic.html

.. _PageRank: http://infolab.stanford.edu/~backrub/google.html

.. _scanpy: https://scanpy.readthedocs.io/en/latest/

.. _scvelo: https://scvelo.readthedocs.io/

.. _documentation: https://cellrank.readthedocs.io

.. _Theislab: https://www.helmholtz-muenchen.de/icb/research/groups/theis-lab/overview/index.html

.. _Peerlab: https://www.mskcc.org/research/ski/labs/dana-pe-er
