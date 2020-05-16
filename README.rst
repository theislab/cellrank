|Travis| |Docs| |Binder| |Codecov|


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

.. |Binder| image:: https://img.shields.io/badge/Tutorial-Pancreas%20Basic-579ACA.svg?logo=data:image/png;base64,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
    :target: https://mybinder.org/v2/gh/theislab/cellrank_notebooks/master?filepath=docs%2Fsource%2Fpancreas_basic.ipynb

.. |Codecov| image:: https://codecov.io/gh/theislab/cellrank/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/theislab/cellrank

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
