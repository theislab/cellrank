.. TODO(michalk8): update the instructions

Installation
============
CellRank requires Python version >= 3.8 to run. We recommend using Miniconda_ to manage the environments.

Anaconda
--------
CellRank can be installed via::

    conda install -c conda-forge -c bioconda cellrank
    # or with extra libraries, useful for large datasets
    conda install -c conda-forge -c bioconda cellrank-krylov

If an error occurs during ``conda install -c conda-forge -c bioconda cellrank-krylov``, please consult the
Dependencies_ section below.

PyPI
----
CellRank is also available on PyPI::

    pip install cellrank

Development Version
-------------------
To stay up-to-date with the newest version, run::

    pip install git+https://github.com/theislab/cellrank@main

Dependencies
------------
Some of the inference tasks that CellRank performs can be broken down to linear algebra problems.
For example, we solve linear systems to compute fate probabilities and we compute partial Schur decompositions and
find metastable states. For these computations to be scalable, we rely on highly optimized libraries which make use
of sparsity structure, parallel implementations and efficient message passing implemented via
`PETSc`_ and `SLEPc`_.
CellRank works without these as well, however, if you would like to apply it to large (>15k cells) datasets,
we recommend you install them.

Below, we give details for installing both `PETSc`_ and `SLEPc`_, in case you've had any issues when installing them
through `PyPI`_ or `Anaconda`_::

    # note: conda alternatives are denoted by alt.

    # update
    sudo apt-get update
    sudo apt-get upgrade

    # install a message passing interface and mpi4py
    sudo apt-get install libopenmpi-dev  # alt.: conda install -c conda-forge openmpi
    pip install --user mpi4py  # alt.: conda install -c anaconda mpi4py

    # install petsc and and petsc4py
    pip install --user petsc  # alt.: conda install -c conda-forge petsc
    pip install --user petsc4py  # alt.: conda install -c conda-forge petsc4py

    # install slepsc and slepsc4py
    pip install --user slepc  # alt.: conda install -c conda-forge slepc
    pip install --user slepc4py  # alt.: conda install -c conda-forge slepc4py

During installation of *petsc*, *petsc4py*, *slepc*, and *slepc4py* the following
error might appear several times::

    ERROR: Failed building wheel for [insert package name here]

but this doesn't matter if the installer finally tells you::

    Successfully installed [insert package name here]

On Mac OS, install `MPICH`_ as a message passing interface and then proceed as above, using either pip or the
installation instructions given on the `PETSc`_ and `SLEPc`_ websites. The `SLEPc`_ homepage even offers a video
tutorial explaining the installation.

If after reading this, you still can't proceed with the installation, feel free to open a `GitHub issue`_.

.. _`Miniconda`: https://conda.pydata.org/miniconda.html
.. _`GitHub issue`: https://github.com/theislab/cellrank/issues/new
.. _`SLEPc`: https://slepc.upv.es/
.. _`PETSc`: https://www.mcs.anl.gov/petsc/
.. _`MPICH`: https://www.mpich.org/
