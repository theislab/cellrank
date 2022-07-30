|PyPI| |Bioconda| |Downloads| |CI| |Notebooks| |Docs| |Codecov|

CellRank: dynamics from multi-view single-cell data
====================================================

.. image:: https://raw.githubusercontent.com/theislab/cellrank/master/resources/images/cellrank_overview.png
   :width: 600px
   :align: center

**CellRank** is a modular framework to study cellular dynamics based on Markov state modeling of multi-view single-cell
data. It estimates differentiation direction based on a varied (and growing!) number of biological priors including RNA
velocity, pseudotime, developmental potential and experimental time points. Read our
:doc:`about cellrank <about_cellrank>` page to learn more.

CellRank scales to large cell numbers, is fully compatible with the `scverse`_ ecosystem, and is easy to use. In the
backend, it is powered by the `pyGPCCA package <https://pygpcca.readthedocs.io/>`_ :cite:`reuter:18,reuter:22`.

CellRank's key applications
----------------------------
- compute initial, terminal and intermediate `macrostates`_ :cite:`reuter:18,reuter:22`.
- infer `fate probabilities`_ towards terminal states.
- visualize `gene expression trends`_ along specific trajectories.
- identify potential `driver genes`_ for each trajectory.
- ... and many more, check out our API.

Getting started with CellRank
------------------------------
We have :doc:`tutorials` and :doc:`examples` that help you getting started; tutorials are longer and explain
computational pipelines,
examples are short and demonstrate individual steps. To learn more about the principles behind CellRank, visit our
:doc:`about CellRank <about_cellrank>` page.

Citing CellRank
----------------
If you find CellRank useful for your research, please visit our :doc:`citing CellRank <citing_cellrank>` page.

Support
--------
We welcome your feedback! Feel free to open an `issue`_ or send
us an `email <mailto:info@cellrank.org>`_ if you encounter a bug, need our help or just want to make a
comment/suggestion.

CellRank in publications
-------------------------
Please check out our manuscript :cite:`lange:22` in Nature Methods to learn more.

Contributing
-------------
We actively encourage any contribution! To get started, please check out the `contribution guide`_.

.. toctree::
    :caption: General
    :maxdepth: 3
    :hidden:

    installation
    api
    developer_api
    about_cellrank
    team
    citing_cellrank
    release_notes
    contributing
    references
    GitHub <https://github.com/theislab/cellrank>
    Discourse <https://discourse.scverse.org/c/ecosytem/cellrank/40>

.. toctree::
    :caption: Gallery
    :maxdepth: 3
    :hidden:

    tutorials
    examples

.. |PyPI| image:: https://img.shields.io/pypi/v/cellrank.svg
    :target: https://pypi.org/project/cellrank
    :alt: PyPI

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/cellrank
    :target: https://bioconda.github.io/recipes/cellrank/README.html
    :alt: Bioconda

.. |Downloads| image:: https://pepy.tech/badge/cellrank
    :target: https://pepy.tech/project/cellrank
    :alt: Downloads

.. |CI| image:: https://img.shields.io/github/workflow/status/theislab/cellrank/Test/master
    :target: https://github.com/theislab/cellrank/actions
    :alt: CI

.. |Notebooks| image:: https://img.shields.io/github/workflow/status/theislab/cellrank_notebooks/CI/master?label=notebooks
    :target: https://github.com/theislab/cellrank_notebooks/actions
    :alt: CI-Notebooks

.. |Docs|  image:: https://img.shields.io/readthedocs/cellrank
    :target: https://cellrank.readthedocs.io/en/stable
    :alt: Documentation

.. |Codecov| image:: https://codecov.io/gh/theislab/cellrank/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/theislab/cellrank
    :alt: Coverage

.. _macrostates: :doc:`notebooks/tutorials/initial_terminal_states`
.. _fate probabilities: :doc:`notebooks/tutorials/fate_probabilities`
.. _driver genes: :doc:`notebooks/tutorials/fate_probabilities`
.. _gene expression trends: :doc:`notebooks/tutorials/gene_trends`
.. _contribution guide: :doc:`contributing`

.. _scverse: https://scverse.org/
.. _Theislab: https://www.helmholtz-muenchen.de/icb/research/groups/theis-lab/overview/index.html
.. _Peerlab: https://www.mskcc.org/research/ski/labs/dana-pe-er
.. _issue: https://github.com/theislab/cellrank/issues/new/choose
