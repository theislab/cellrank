|PyPI| |Bioconda| |Downloads| |Discourse| |CI| |Notebooks| |Docs| |Codecov|

CellRank for directed single-cell fate mapping
==============================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/main/docs/_static/img/cellrank_overview.png
    :width: 600px
    :align: center

.. sidebar:: Key Contributors

    * `Marius Lange`_: lead developer, initial CellRank conception, maintainer
    * `Michal Klein`_: senior developer, design & architecture, maintainer
    * `Philipp Weiler`_: developer

.. _Marius Lange: https://twitter.com/MariusLange8
.. _Michal Klein: https://github.com/michalk8
.. _Philipp Weiler: https://twitter.com/PhilippWeiler7

**CellRank** is a toolkit to uncover cellular dynamics based on Markov state modeling of single-cell data. It contains
two main modules: `kernels`_ compute cell-cell transition probabilities and `estimators`_ generate hypothesis based on
these. Our kernels work with a variety of input data including `RNA velocity`_ (see `La Manno et al. (2018)`_ and
`Bergen et al. (2020)`_), `cellular similarity`_ (both transcriptomic and spatial) and `pseudotime`_, among others.
Our `VelocityKernel`_ takes into account **uncertainty in the velocities** and allows you to aggregate the short-range
fate relations given by RNA velocity into longer trends along the phenotypic manifold. Our main estimator is
*Generalized Perron Cluster Cluster Analysis* (`Reuter et al. (2018)`_) which coarse-grains the Markov chain
into a set of macrostates which represent initial, terminal and intermediate states. For each transient cell,
we compute its fate probability towards any terminal state. We show an example of such a fate map in the figure above,
which has been computed using the data of `pancreatic endocrinogenesis`_. CellRank combines `kernels`_ and `estimators`_
with a powerful `plotting API`_, enabling you to visualize e.g. smooth `gene expression trends`_ along lineages or
fate-informed `circular embeddings`_, to name just a few.

CellRank scales to large cell numbers, is fully compatible with `scanpy`_ and `scvelo`_ and is easy to use.
For **installation instructions**, **documentation** and **tutorials**, visit `cellrank.org`_.

Manuscript
^^^^^^^^^^
Please check out our manuscript `Lange et al. (2022)`_ in **Nature Methods** to learn more.

Getting started with CellRank
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you're new to CellRank, make sure to go though the `basic tutorial`_ which introduces you to CellRank's high-level
API. Most biological systems require a bit more control, so be sure to check out the `kernels and estimators tutorial`_
which allows to unlock the full power of CellRank. If you want to see individual functions in action, visit our
`gallery`_.

CellRank's key applications
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- compute initial & terminal as well as intermediate `macrostates`_ of your biological system
- infer `fate probabilities`_ towards the terminal states for each individual cell
- visualize `gene expression trends`_ along specific lineages while accounting for the continuous nature of
  fate determination
- identify potential `driver genes`_ for each identified cellular trajectory

Installation
^^^^^^^^^^^^
Install CellRank by running::

    conda install -c conda-forge -c bioconda cellrank
    # or with extra libraries, useful for large datasets
    conda install -c conda-forge -c bioconda cellrank-krylov

or via PyPI::

    pip install cellrank

Why is it called "CellRank"?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CellRank **does not** rank cells, we gave the package this name because just like Google's original `PageRank`_
algorithm, it works with Markov chains to aggregate relationships between individual objects (cells vs. websites)
to learn about more global properties of the underlying dynamics (initial & terminal states and fate probabilities vs.
website relevance).

Support
^^^^^^^
We welcome your feedback! Feel free to open an `issue <https://github.com/theislab/cellrank/issues/new/choose>`_, send
us an `email <mailto:info@cellrank.org>`_ or `tweet`_ if you encounter a bug, need our help or just want to make a
comment/suggestion.

Contributing
^^^^^^^^^^^^
We actively encourage any contribution! To get started, please check out both the `contribution guide`_ as well as the
`external API`_. CellRank's modular structure makes it easy to contribute, be it a new method to compute cell-cell
transition probabilities (`kernels`_), a new way to analyze a transition matrix (`estimators`_) or an addition to the
`plotting API`_. If you're thinking of contributing a new kernel, we have a `kernel tutorial`_ that guides you trough
the process.

CellRank was developed in collaboration between the `Theislab`_ and the `Peerlab`_.

.. |PyPI| image:: https://img.shields.io/pypi/v/cellrank.svg
    :target: https://pypi.org/project/cellrank
    :alt: PyPI

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/cellrank
    :target: https://anaconda.org/bioconda/cellrank
    :alt: Bioconda

.. |Downloads| image:: https://pepy.tech/badge/cellrank
    :target: https://pepy.tech/project/cellrank
    :alt: Downloads

.. |Discourse| image:: https://img.shields.io/discourse/posts?color=yellow&logo=discourse&server=https%3A%2F%2Fdiscourse.scverse.org
    :target: https://discourse.scverse.org/
    :alt: Discourse

.. |CI| image:: https://img.shields.io/github/actions/workflow/status/theislab/cellrank/test.yml?branch=main
    :target: https://github.com/theislab/cellrank/actions
    :alt: CI

.. |Notebooks| image:: https://img.shields.io/github/actions/workflow/status/theislab/cellrank_notebooks/ci.yml?branch=main&label=notebooks
    :target: https://github.com/theislab/cellrank_notebooks/actions
    :alt: CI-Notebooks

.. |Docs|  image:: https://img.shields.io/readthedocs/cellrank
    :target: https://cellrank.readthedocs.io/en/stable
    :alt: Documentation

.. |Codecov| image:: https://codecov.io/gh/theislab/cellrank/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/theislab/cellrank
    :alt: Coverage


.. _La Manno et al. (2018): https://doi.org/10.1038/s41586-018-0414-6
.. _Bergen et al. (2020): https://doi.org/10.1038/s41587-020-0591-3
.. _Reuter et al. (2018): https://doi.org/10.1021/acs.jctc.8b00079
.. _pancreatic endocrinogenesis: https://doi.org/10.1242/dev.173849
.. _cellrank.org: https://cellrank.org

.. _kernels: https://cellrank.readthedocs.io/en/stable/classes.html#kernels
.. _estimators: https://cellrank.readthedocs.io/en/stable/classes.html#estimators
.. _plotting API: https://cellrank.readthedocs.io/en/stable/api.html#module-cellrank.pl
.. _external API: https://cellrank.readthedocs.io/en/stable/external_api.html
.. _contribution guide: https://github.com/theislab/cellrank/blob/main/CONTRIBUTING.rst

.. _RNA velocity: https://cellrank.readthedocs.io/en/stable/classes.html#velocity-kernel
.. _VelocityKernel: https://cellrank.readthedocs.io/en/stable/classes.html#velocity-kernel
.. _cellular similarity: https://cellrank.readthedocs.io/en/stable/classes.html#connectivity-kernel
.. _pseudotime: https://cellrank.readthedocs.io/en/stable/classes.html#pseudotime-kernel

.. _gene expression trends: https://cellrank.readthedocs.io/en/stable/api/cellrank.pl.gene_trends.html#cellrank.pl.gene_trends
.. _circular embeddings: https://cellrank.readthedocs.io/en/stable/api/cellrank.pl.circular_projection.html

.. _basic tutorial: https://cellrank.readthedocs.io/en/stable/cellrank_basics.html
.. _kernel tutorial: https://cellrank.readthedocs.io/en/stable/creating_new_kernel.html
.. _kernels and estimators tutorial: https://cellrank.readthedocs.io/en/stable/kernels_and_estimators.html

.. _scanpy: https://scanpy.readthedocs.io/en/stable/
.. _scvelo: https://scvelo.readthedocs.io/

.. _Theislab: https://www.helmholtz-muenchen.de/icb/research/groups/theis-lab/overview/index.html
.. _Peerlab: https://www.mskcc.org/research/ski/labs/dana-pe-er
.. _`tweet`: https://twitter.com/MariusLange8
.. _Lange et al. (2022): https://www.nature.com/articles/s41592-021-01346-6
.. _PageRank: https://en.wikipedia.org/wiki/PageRank#cite_note-1

.. _gallery: https://cellrank.readthedocs.io/en/stable/auto_examples/index.html
.. _macrostates: https://cellrank.readthedocs.io/en/stable/auto_examples/estimators/compute_macrostates.html
.. _fate probabilities: https://cellrank.readthedocs.io/en/stable/auto_examples/estimators/compute_abs_probs.html
.. _driver genes: https://cellrank.readthedocs.io/en/stable/auto_examples/estimators/compute_lineage_drivers.html
