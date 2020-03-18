.. Cellrank documentation master file, created by
   sphinx-quickstart on Sun Dec  1 21:13:53 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

|Travis| |License|


CellRank - Probabilistic Trajectory Inference based on RNA Velocity
===================================================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/index_figure.png
   :width: 600px
   :align: center

**CellRank** is a toolkit to uncover cellular development based on scRNA-seq
data with RNA Velocity annotation, see [Manno18]_ and [Bergen19]_.

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

Report issues and see the code on `GitHub <https://github.com/theislab/cellrank>`__.

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

   Pancreas <https://cellrank-notebooks.readthedocs.io/en/latest/pancreas_basic.html>


.. |Travis| image:: https://travis-ci.org/theislab/cellrank.svg?branch=master
    :target: https://travis-ci.org/theislab/cellrank

.. |License| image:: https://img.shields.io/github/license/theislab/cellrank
