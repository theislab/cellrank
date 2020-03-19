|Travis| |Docs| |License|


CellRank - Probabilistic Trajectory Inference based on RNA Velocity
===================================================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/index_figure.png
   :width: 600px
   :align: center

**CellRank** is a toolkit to uncover cellular development based on scRNA-seq
data with RNA Velocity annotation, see `La Manno et al. (2018) <https://doi.org/10.1038/s41586-018-0414-6>`_
and `Bergen et al. (2019) <https://doi.org/10.1101/820936>`_.

CellRank utilises the time derivative of gene expression given by RNA Velocity
to construct a Markov Chain. The information given by RNA Velocity is combined
with transcriptomic similarity and density corrected to yield a robust estimate
of cellular development directly in high dimensional gene expression space.

Based on this Markov Chain, we infer start and endpoints of development as well
as lineage probabilities, i.e. for each cell, the probability of it reaching
any of the inferred endpoints of development. CellRank offers many possibilities
to utilise these lineage probabilities in downstream analysis, e.g. to investigate
the fate of early cell clusters or to plot gene expression trends along a given
lineage.

CellRank is fully compatible with `scanpy <https://scanpy.readthedocs.io/en/latest/>`_
and `scvelo <https://scvelo.readthedocs.io/>`_.

CellRank requires Python3 version >= 3.6 and can be installed via::

    pip install git+https://github.com/theislab/cellrank

See the documentation at `<https://cellrank.readthedocs.io>`_, which
includes tutorial notebooks and a description of our complete API.


.. |Travis| image:: https://travis-ci.org/theislab/cellrank.svg?branch=master
    :target: https://travis-ci.org/theislab/cellrank

.. |Docs|  image:: https://img.shields.io/readthedocs/cellrank

.. |License| image:: https://img.shields.io/github/license/theislab/cellrank
