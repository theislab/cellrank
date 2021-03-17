Contributing guide
==================

Table of Contents
=================
- `Introduction`_

  - `Setup`_
  - `Codebase structure`_

- `Adding new features`_

  - `Creating new kernel`_
  - `Creating new estimator`_
  - `Creating new model`_
  - `Adding external package`_

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
The Squidpy project:

- `cellrank <cellrank>`_: the root of the package.

  - `cellrank/datasets <cellrank/datasets>`_: contains downloadable datasets.
  - `cellrank/external <cellrank/external>`_: contains all external kernels/estimators/models.
  - `cellrank/pl <cellrank/pl>`_: contains high level plotting functions.
  - `cellrank/tl/kernels <cellrank/tl/kernels>`_: contains kernels that create transition matrices.
  - `cellrank/tl/estimators <cellrank/tl/estimators>`_: contains estimators which estimate fate probabilities based on kernels.
  - `cellrank/ul/models <cellrank/ul/models>`_: contains smooth gene expression models.

Adding new features
~~~~~~~~~~~~~~~~~~~

Creating new kernel
-------------------

Creating new estimator
----------------------

Creating new model
------------------

Adding external package
-----------------------

Running tests
~~~~~~~~~~~~~
We use ``tox`` alongside with ``pytest`` to automate the testing process. To run the tests, use the command below::

    tox -e py{37,38,39}-{linux,macos}

depending on the Python version(s) in your ``PATH`` and your operating system. For example, on Linux and Python3.8,
you would run::

    tox -e py38-linux

To run only a subset of tests, run::

    tox -e <environment> -- <name>

where ``<name>`` can be a path to a test file/directory or a name of a test function/class. In the example below, we'd
run only tests in the file

    tox -e py38-linux -- tests/graph/test_kernels.py

Building documentation
~~~~~~~~~~~~~~~~~~~~~~
In order to build the documentation, run one of the commands below,
depending on whether you also want to build the examples::

    tox -e docs  # builds examples as well, takes longer
    tox -e shallow-docs  # does not build the examples

If you need to clean the artifacts from previous documentation builds, run::

    tox -e clean-docs

Troubleshooting
---------------
- **problems with running tox commands**
  Try recreating the environment as::

    tox -e <environment> --recreate

  If this didn't work, you can purge the whole ``.tox`` directory as ``rm -rf .tox``.

- **can't commit because of pre-commit**
  Sometimes, the quirks of the linter ... You can temporarily bypass pre-commit by commiting as::

    git commit --no-verify
