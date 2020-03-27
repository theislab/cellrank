|Travis| |Docs| |Codecov| |License|


CellRank - Continuous Lineage Decisions Uncovered by RNA Velocity
===================================================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/index_figure_endpoints.png
   :width: 600px
   :align: center

**CellRank** is a toolkit to uncover cellular development based on scRNA-seq data with RNA velocity annotation,
see [Manno18]_ and [Bergen19]_.
|br|
In the figure above, we show the main features of CellRank applied to [Panc19]_ - starting from
RNA velocities **(a)**, we infer root cells **(b)** and final cells **(c)**, and we compute how likely each cell
is to develop towards each of the identified groups of final cells **(d)**.
|br|
See our `tutorial`_ to learn how to apply CellRank to your own data.

CellRank utilizes the time derivative of gene expression given by RNA velocity
to construct a Markov chain. The information given by RNA velocity is combined
with transcriptomic similarity and density corrected to yield a robust estimate
of cellular development directly in high dimensional gene expression space.
|br|
CellRank obtained its name due to conceptual similarities with `PageRank`_.
Google's algorithm for ranking web pages. Both algorithms construct a Markov chain
and use spectral methods to study the long term evolution of this process.

Based on this Markov chain, we infer root and final cells of development as well
as lineage probabilities, i.e. for each cell, the probability of it reaching
any of the inferred root and final cells of development. CellRank offers many possibilities
to utilize these lineage probabilities in downstream analysis, e.g. to investigate
the fate of early cell clusters or to plot gene expression trends along a given lineage.

CellRank is fully compatible with `scanpy`_. and `scvelo`_
and was developed as a collaboration between `Theislab`_ and `Peerlab`_.


.. toctree::
    :caption: General
    :maxdepth: 2
    :hidden:

    installation
    api
    classes
    references

.. toctree::
   :caption: Tutorials
   :maxdepth: 2
   :hidden:

   Pancreas Basic <https://cellrank-notebooks.readthedocs.io/en/latest/pancreas_basic.html>
   Pancreas Advanced <https://cellrank-notebooks.readthedocs.io/en/latest/pancreas_advanced.html>


.. |Travis| image:: https://travis-ci.org/theislab/cellrank.svg?branch=master
    :target: https://travis-ci.org/theislab/cellrank

.. |Docs|  image:: https://img.shields.io/readthedocs/cellrank

.. |Codecov| image:: https://codecov.io/gh/theislab/cellrank/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/theislab/cellrank

.. |License| image:: https://img.shields.io/github/license/theislab/cellrank

.. _tutorial: <https://cellrank-notebooks.readthedocs.io/en/latest/pancreas_basic.html>:

.. _PageRank: <http://infolab.stanford.edu/~backrub/google.html>:

.. _documentation: <https://cellrank.readthedocs.io>:

.. _scanpy: <https://scanpy.readthedocs.io/en/latest/>:

.. _scvelo: <https://scvelo.readthedocs.io/>:

.. _Theislab: <https://www.helmholtz-muenchen.de/icb/research/groups/theis-lab/overview/index.html>:

.. _Peerlab: <https://www.mskcc.org/research/ski/labs/dana-pe-er>:

.. |br| raw:: html

  <br/>
