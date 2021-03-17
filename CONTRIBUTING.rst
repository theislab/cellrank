Contributing guide
==================

Table of Contents
=================
- `Introduction`_

  - `Setup`_
  - `Codebase structure`_

- `Adding new features`_

  - `Adding external tools`_
  - `Creating new kernel`_

- `Running tests`_
- `Building documentation`_
- `Troubleshooting`_

Introduction
~~~~~~~~~~~~

Setup
-----
Clone CellRank from source as::

    git clone https://github.com/theislab/cellrank
    cd cellrank

Install development requirements in editable mode, as well as pre-commit::

    pip install -e'.[dev]'
    pre-commit install

Codebase structure
------------------
The CellRank is structured as follows:

- `cellrank <cellrank>`_: the root of the package.

  - `cellrank/datasets <cellrank/datasets>`_: contains datasets available for download.
  - `cellrank/external <cellrank/external>`_: contains all external kernels/estimators/models.
  - `cellrank/pl <cellrank/pl>`_: contains high level plotting functions.
  - `cellrank/tl/kernels <cellrank/tl/kernels>`_: contains classes that create transition matrices.
  - `cellrank/tl/estimators <cellrank/tl/estimators>`_: contains classes which estimate fate probabilities based on kernels.
  - `cellrank/ul/models <cellrank/ul/models>`_: contains classes that compute smoothed gene expression trends.

- `tests <tests>`_: unit tests, see `Running tests`_ for more information.
- `docs <docs>`_: documentation, see `Building documentation`_ for more.
- `examples <examples>`_: the examples as seen in the
  `documentation <https://cellrank.readthedocs.io/en/latest/auto_examples/index.html>`_.

Adding new features
~~~~~~~~~~~~~~~~~~~
We welcome every contribution to CellRank, whether it's a new kernel/estimator/model,
or a new high-level function/modified functionality. Below, we mostly focus on kernels.

Creating new kernel
-------------------
If you'd like to contribute a new kernel to CellRank, either as a part of an external package or not,
we have created a short, 5 minute `tutorial <https://cellrank.readthedocs.io/en/latest/creating_new_kernel.html>`_
that shows you how to do it.

Adding external tools
---------------------
If you already have/know of an external package that could be made available through CellRank, we'd love to hear!

In short, external packages should be accessible through ``cellrank.external`` and a thin wrapper code should be placed
under:

- `cellrank/external/kernels <cellrank/external/kernels>`_: if adding an external kernel. For example,
  see `OTKernel <cellrank/external/kernels/_statot_kernel.py>`_ from `statOT <https://github.com/zsteve/StationaryOT>`_.
- `cellrank/external/estimators <cellrank/external/estimators>`_: if adding an external estimator.
- `cellrank/external/models <cellrank/external/models>`_: if adding an external model.

You can either open an `issue <https://github.com/theislab/cellrank/issues/new/choose>`_ with a suggestion or
directly submit a `PR <https://github.com/theislab/cellrank/pulls>`_ containing the new addition.


Running tests
~~~~~~~~~~~~~
We use ``tox`` alongside with ``pytest`` to automate the testing process. To run the tests, use the command below::

    tox -e py{37,38,39}-{linux,macos}

depending on the Python version(s) in your ``PATH`` and your operating system. For example, on Linux and Python 3.8,
you would run::

    tox -e py38-linux

To run only a subset of tests, run::

    tox -e <environment> -- <name>

where ``<name>`` can be a path to a test file/directory or a name of a test function/class. In the example below, we'd
run only tests in the file::

    tox -e py38-linux -- tests/test_kernels.py

Building documentation
~~~~~~~~~~~~~~~~~~~~~~
In order to build the documentation, run one of the commands below,
depending on whether you also want to build the examples::

    tox -e docs  # builds examples as well, takes longer
    tox -e shallow-docs  # does not build the examples

If you need to clean the artifacts from previous documentation builds, run::

    tox -e clean-docs

Troubleshooting
~~~~~~~~~~~~~~~
- **I have problems with running some tox commands**

  Try recreating the environment as::

    tox -e <environment> --recreate

  If this didn't work, you can purge the whole ``.tox`` directory as ``rm -rf .tox``.

- **I can't commit because of pre-commit**

  Sometimes, it can be hard to satisfy the linting step. You can temporarily bypass it by committing as::

    git commit --no-verify
