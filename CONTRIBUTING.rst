Contributing guide
==================

Table of Contents
=================
- `Introduction`_

  - `Setup`_
  - `Codebase structure`_
  - `Branching structure`_

- `Adding new features`_

  - `Adding external tools`_
  - `Creating a new kernel`_

- `Running tests`_
- `Documentation`_

  - `Building documentation`_
  - `Writing documentation`_

- `Tutorials/examples`_
- `Making use of GitHub issues/discussions`_
- `Maintainer notes`_

  - `Making a new release`_

- `Troubleshooting`_

Introduction
~~~~~~~~~~~~

Setup
-----
Clone CellRank from source as::

    git clone https://github.com/theislab/cellrank
    cd cellrank

Install development requirements in editable mode, as well as pre-commit as::

    pip install -e'.[dev]'
    pre-commit install

Codebase structure
------------------
The CellRank is structured as follows:

- `cellrank <cellrank>`_: the root of the package.

  - `cellrank/datasets <cellrank/datasets>`__: contains datasets that are available for download.
  - `cellrank/external <cellrank/external>`_: contains all external kernels/estimators/models.
  - `cellrank/pl <cellrank/pl>`_: the plotting module, containing high level plotting functions.
  - `cellrank/tl <cellrank/tl>`_: the tools module, containing high-level tool functions, linear solvers, etc.

    - `cellrank/tl/kernels <cellrank/tl/kernels>`_: contains classes that create transition matrices based on
      various input features, such as transcriptomic or spatial similarities, RNA velocity, pseudotime, etc.
    - `cellrank/tl/estimators <cellrank/tl/estimators>`_: contains classes with perform computations on Markov chains
      defined through the transition matrices obtained in the previous step through one or several kernels.
      These computations usually involve high-level abstractions of the Markov chains, such as coarse-graining
      into macrostates, finding initial and terminal states, estimating fate probabilities, etc.

  - `cellrank/ul <cellrank/ul>`_: the utilities modules, containing mostly models for gene trend smoothing,
    parallelization, etc.

    - `cellrank/ul/models <cellrank/ul/models>`_: contains classes that compute smoothed gene expression trends.

- `tests <tests>`_: unit tests, see `Running tests`_ for more information.
- `docs <docs>`_: documentation, see `Documentation`_ for more.
- `examples <examples>`__: the examples as seen in the
  `documentation <https://cellrank.readthedocs.io/en/latest/auto_examples/index.html>`__, see `Tutorials/examples`_
  for more information.

Branching structure
~~~~~~~~~~~~~~~~~~~
We use sightly modified branching structured described
`here <https://nvie.com/posts/a-successful-git-branching-model/>`_, in short:

- ``feature/X``: for features, enhancements, etc.; it is based off the ``dev`` branch and squash merged into it again.
- ``fix/X``: for bug fixes; is based off the ``dev`` branch and squash merged into it again.
- ``release/vX.X.X``: prior to making a release; it is based off ``dev`` and then merged (using merge commit) into
  ``master`` and ``dev``, after bumping the version and making release notes.
  Both are usually done through CI, see `Making a new release`_.

Adding new features
~~~~~~~~~~~~~~~~~~~
We welcome every contribution to CellRank, whether it's a new kernel/estimator/model,
or a new high-level function/modified functionality. Below, we mostly focus on kernels.

Creating a new kernel
---------------------
Creating a new kernel is the primary way to interface with CellRank if you have created a method which computes a
transition matrix among single-cells. This can be based on any input features that are measured in an experiment, such
as transcriptional similarity, chromatin accessibility, RNA velocity, experimental timepoints, metabolic labeling, etc.

If you have already implemented your method in a separate package, it can be wrapped in a CellRank kernel which will
prompt the user to install you package upon initialization, hence boosting your package download numbers.
Making your package available through a CellRank kernel will give it more visibility and allow users to benefit
from all of CellRank's downstream functionalities available through estimators.

If you'd like to contribute a new kernel to CellRank, either as a part of an external package or not,
we have created a short, 5 minute `tutorial <https://cellrank.readthedocs.io/en/latest/creating_new_kernel.html>`_
that shows you how to do it.

Adding external tools
---------------------
If you already have/know of an external package that could be made available through CellRank, we'd love to hear!

In short, external packages should be accessible through ``cellrank.external`` and a thin wrapper code should be placed
under:

- `cellrank/external/kernels <cellrank/external/kernels>`_: if adding an external kernel. As an example, see the
  `OTKernel <cellrank/external/kernels/_statot_kernel.py>`_ from `statOT <https://github.com/zsteve/StationaryOT>`_.
- `cellrank/external/estimators <cellrank/external/estimators>`_: if adding an external estimator.
- `cellrank/external/models <cellrank/external/models>`_: if adding an external model.

You can either open an `issue <https://github.com/theislab/cellrank/issues/new/choose>`__ with a suggestion or
directly submit a `PR <https://github.com/theislab/cellrank/pulls>`_ containing the addition.
Furthermore, external package dependencies should be added to `setup.py's <setup.py>`_ ``'external'``
in ``extras_require`` and a corresponding entry should be made in the external API
`documentation <docs/source/external_api.rst>`__.

Running tests
~~~~~~~~~~~~~
We use ``tox`` alongside with ``pytest`` to automate the testing process. To run the tests, use the command below::

    tox -e py{37,38,39}-{linux,macos}

depending on the Python version(s) in your ``PATH`` and your operating system. For example, on Linux and Python 3.8,
you would run::

    tox -e py38-linux

Note that during the first invocation, it can take several minutes to start the testing. This is due to the fact that
PETSc/SLEPc needs to be built. To run only a subset of tests, run::

    tox -e <environment> -- -k <name>

where ``<name>`` can be a path to a test file/directory or a name of a test function/class, e.g.::

    tox -e py38-linux -- tests/test_kernels.py  # run tests in the speciied file
    tox -e py38-linux -- -k "TestGPCCA"  # run tests grouped in the `TestGPCCA` class

Documentation
~~~~~~~~~~~~~

Building documentation
----------------------
In order to build the documentation, run one of the commands below, depending on whether you also want to build the
examples::

    tox -e docs  # builds the examples, takes longer (~10 mins)
    tox -e shallow-docs  # does not build the examples

If you need to clean the artifacts from previous documentation builds, run::

    tox -e clean-docs

Writing documentation
---------------------
We use ``numpy``-style docstrings for the documentation with the following additions and modifications:

- no type hints in the docstring (optionally applies also for the return statement) should be used,
  since all functions are required to have the type hints in their signatures.
- when referring to some argument within the same docstring, enclose that reference in \`\`.
- when referring to an argument of a class from within that class, use ``:paramref:`attribute```.
- optional, but recommended: when referring to attributes of a foreign class, use ``:attr:`qualified_name```, such as
  ``:attr:`anndata.AnnData.obs```.
- use ``docrep`` for repeating documentation.

Below is an example of how a docstring should look like::

    from cellrank.ul._docs import d

    @d.dedent  # using docrep to interpolate %(adata)s
    def some_function(adata: AnnData, key: str) -> float:
        """
        This is a short one-line header.

        Here you can add multi-paragraph explanation, if needed.

        Parameters
        ----------
        %(adata)s
        key
            Some key in :attr:`anndata.AnnData.obs`.

        Returns
        --------
        Some return description.
        """


Making use of GitHub issues/discussions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Opening an `issue <https://github.com/theislab/cellrank/issues>`__ or
starting a `discussion <https://github.com/theislab/cellrank/discussions>`_ is the primary way to get help.
Issues are used mostly for feature requests or for fixing bugs, whereas discussions can be used to ask conceptual
questions, algorithmic/biological questions or just to exchange ideas.

Maintainer notes
~~~~~~~~~~~~~~~~

Making a new release
--------------------
New release is always created when a new tag is pushed to GitHub. When that happens, a new CI job starts the
testing machinery. If all the tests pass, new release will be created on PyPI. Bioconda will automatically notice that
a new release has been made and an automatic PR will be made to
`bioconda-recipes <https://github.com/bioconda/bioconda-recipes/pulls>`_.
Extra care has to be taken when updating runtime dependencies - this is not automatically picked up by Bioconda
and a separate PR with the updated ``recipe.yaml`` will have to be made.

Easiest way to create a new release it to create a branch named ``release/vX.X.X`` and push it onto GitHub. The CI
will take care of the following:

- create the new release notes (and empty the ``dev`` release notes)
- bump the version and create a new tag
- run tests on the ``release/vX.X.X`` branch
- publish on PyPI after the tests have passed
- merge ``release/vX.X.X`` into ``master`` and ``dev``

Alternatively, it's possible to create a new release using ``bump2version``, which can be installed as::

    pip install bump2version

Depending on what part of the version you want to update, you can run on ``master``::

    bump2version {major,minor,patch}

By default, this will create a new tagged commit, automatically update the ``__version__`` wherever necessary.
Afterwards, you can just push the changes to upstream by running::

    git push --atomic <branch> <tag>

or set ``push.followtags=true`` in your git config and do a regular ``git push``. In this case, CI will not run
create any release notes, run tests or do any merges.

Tutorials/examples
~~~~~~~~~~~~~~~~~~
While our tutorials focus on an entire workflow or module of CellRank, i.e. using RNA velocity and similarity
to compute terminal states, examples focus on a single function/method and show how it can be used in practice.

The tutorials are hosted in a separate `repo <https://github.com/theislab/cellrank_notebooks>`_, whereas examples
are hosted in this repo, under `examples <examples>`__. Both tutorials and examples use already preprocessed datasets
from `cellrank/datasets <cellrank/datasets>`__, with precomputed attributes, such as velocities, pseudotime, etc.

If you wish to contribute your own example (e.g. for an external kernel), you just need to write a ``.py`` file, similar
to `this one <examples/other/compute_kernel_tricks.py>`_.
The filenames should be prefixed with either ``compute_`` or ``plot_``, depending on what they do, i.e. whether they
show a computational or a plotting functionality.

Tutorials, on the other hand, are written as Jupyter notebooks. However, they are still tested on the CI to make sure
they run properly with the newest version of CellRank. Since they require more effort to create than the examples,
it's best to first start a new issue/discussion before adding them, see also `Making use of GitHub issues/discussions`_.

Troubleshooting
~~~~~~~~~~~~~~~
- **I have problems with running some tox commands**

  Try recreating the environment as::

    tox -e <environment> --recreate

  If this didn't work, you can purge the whole ``.tox`` directory as ``rm -rf .tox``.

- **I can't commit because of pre-commit**

  Sometimes, it can be hard to satisfy the linting step. You can temporarily bypass it by committing as::

    git commit --no-verify

- **I have an issue which this section does not solve**

  Please see `Making use of GitHub issues/discussions`_ on how to create a new issue or how to start a discussion.
