|PyPI| |Bioconda| |Downloads| |CI| |Notebooks| |Docs| |Codecov|

CellRank: dynamics from multi-view single-cell data
===================================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/main/docs/source/_static/img/cellrank_overview.png
    :width: 600px
    :align: center

**CellRank** :cite:`lange:22` is a modular framework to study cellular dynamics based on Markov state modeling of
multi-view single-cell data. It estimates differentiation direction based on a varied (and growing!) number of
biological priors including RNA velocity, pseudotime, developmental potential and experimental time points.
See :doc:`about CellRank <about/index>` to learn more.

CellRank scales to large cell numbers, is fully compatible with the `scverse`_ ecosystem, and is easy to use. In the
backend, it is powered by the `pyGPCCA package <https://github.com/msmdev/pyGPCCA>`_ :cite:`reuter:19,reuter:22`. Feel
free to open an `issue`_ or send us an `email <mailto:info@cellrank.org>`_ if you encounter a bug, need our help or just
want to make a comment/suggestion.

CellRank's key applications
---------------------------
- Under reconstruction. .. TODO(Marius1311): please re-factor this page

Getting started with CellRank
-----------------------------
We have :doc:`notebooks/tutorials/index` that will help you getting started. To see CellRank in action, explore our
manuscript :cite:`lange:22` in Nature Methods.

Contributing
------------
We actively encourage any contribution! To get started, please check out the :doc:`contributing`.

.. toctree::
    :caption: General
    :maxdepth: 3
    :hidden:

    installation
    api/index
    notebooks/tutorials/index
    release_notes
    contributing
    references

.. toctree::
    :caption: About
    :maxdepth: 3
    :hidden:

    about/index
    about/cite
    about/team
    GitHub <https://github.com/theislab/cellrank>
    Discourse <https://discourse.scverse.org/c/ecosytem/cellrank/40>

.. |PyPI| image:: https://img.shields.io/pypi/v/cellrank.svg
    :target: https://pypi.org/project/cellrank
    :alt: PyPI

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/cellrank
    :target: https://anaconda.org/bioconda/cellrank
    :alt: Bioconda

.. |Downloads| image:: https://pepy.tech/badge/cellrank
    :target: https://pepy.tech/project/cellrank
    :alt: Downloads

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

.. _scverse: https://scverse.org/
.. _issue: https://github.com/theislab/cellrank/issues/new/choose
