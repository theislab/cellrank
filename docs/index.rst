|PyPI| |Downloads| |CI| |Docs| |Codecov| |Discourse|

CellRank 2: Unified fate mapping in multiview single-cell data
==============================================================
.. image:: _static/img/light_mode_overview.png
    :width: 600px
    :align: center
    :class: only-light

.. image:: _static/img/dark_mode_overview.png
    :width: 600px
    :align: center
    :class: only-dark

**CellRank** :cite:`lange:22,weiler:23` is a modular framework to study cellular dynamics based on Markov state modeling of
multi-view single-cell data. See :doc:`about CellRank <about/index>` to learn more and :doc:`here <about/cite>` for how to correctly
cite our work.

CellRank scales to large cell numbers, is fully compatible with the `scverse`_ ecosystem, and is easy to use. In the
backend, it is powered by the `pyGPCCA package <https://github.com/msmdev/pyGPCCA>`_ :cite:`reuter:19,reuter:22`. Feel
free to open an `issue`_ or send us an `email <mailto:info@cellrank.org>`_ if you encounter a bug, need our help or
just want to make a comment/suggestion.

.. important::
    If you're moving from CellRank 1 to CellRank 2, check out :doc:`../about/version2`.

CellRank's Key Applications
---------------------------
- Estimate differentiation direction based on a varied number of biological priors, including
  :doc:`pseudotime <notebooks/tutorials/kernels/300_pseudotime>`,
  :doc:`developmental potential <notebooks/tutorials/kernels/400_cytotrace>`,
  :doc:`RNA velocity <notebooks/tutorials/kernels/200_rna_velocity>`,
  :doc:`experimental time points <notebooks/tutorials/kernels/500_real_time>`, and :mod:`more <cellrank.kernels>`.
- Compute initial, terminal and intermediate :doc:`macrostates <notebooks/tutorials/estimators/600_initial_terminal>`
  :cite:`reuter:19,reuter:22`.
- Infer :doc:`fate probabilities and driver genes <notebooks/tutorials/estimators/700_fate_probabilities>`.
- Visualize and cluster :doc:`gene expression trends <notebooks/tutorials/estimators/800_gene_trends>`.
- ... and much more, check out our :doc:`API <api/index>`.

Getting Started with CellRank
-----------------------------
We have :doc:`notebooks/tutorials/index` to help you getting started. To see CellRank in action, explore our
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
    about/version2
    about/cite
    about/team
    GitHub <https://github.com/theislab/cellrank>
    Discourse <https://discourse.scverse.org/c/ecosytem/cellrank/40>

.. |PyPI| image:: https://img.shields.io/pypi/v/cellrank.svg
    :target: https://pypi.org/project/cellrank
    :alt: PyPI

.. |Downloads| image:: https://pepy.tech/badge/cellrank
    :target: https://pepy.tech/project/cellrank
    :alt: Downloads

.. |Discourse| image:: https://img.shields.io/discourse/posts?color=yellow&logo=discourse&server=https%3A%2F%2Fdiscourse.scverse.org
    :target: https://discourse.scverse.org/c/ecosystem/cellrank/
    :alt: Discourse

.. |CI| image:: https://img.shields.io/github/actions/workflow/status/theislab/cellrank/test.yml?branch=main
    :target: https://github.com/theislab/cellrank/actions
    :alt: CI

.. |Docs|  image:: https://img.shields.io/readthedocs/cellrank
    :target: https://cellrank.readthedocs.io/
    :alt: Documentation

.. |Codecov| image:: https://codecov.io/gh/theislab/cellrank/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/theislab/cellrank
    :alt: Coverage

.. _scverse: https://scverse.org/
.. _issue: https://github.com/theislab/cellrank/issues/new/choose
