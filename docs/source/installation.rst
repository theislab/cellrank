Installation
============
CellRank requires Python version >= 3.6 to run. We recommend using Miniconda_ to manage the environments.

PyPI
~~~~
Install CellRank by running::

    pip install cellrank
    # or with highly optimized libraries - this can take a long time
    pip install cellrank[krylov]

If an error occurs during ``pip install cellrank[krylov]``, please consult the Dependencies_ section below.

Anaconda
~~~~~~~~
CellRank is also available on Anaconda and can be installed via::

    conda install -c conda-forge -c bioconda cellrank
    # or with highly optimized libraries - recommended approach
    conda install -c conda-forge -c bioconda cellrank-krylov

If an error occurs during ``conda install -c conda-forge -c bioconda cellrank-krylov``, please consult the
Dependencies_ section below.

Development Version
~~~~~~~~~~~~~~~~~~~
To stay up-to-date with the newest version, run::

    git clone https://github.com/theislab/cellrank
    cd cellrank
    pip install -e '.[dev]'
    python-vendorize

``-e`` stands for ``--editable`` and makes sure that your environment is updated
when you pull new changes from GitHub. The ``'[dev]'`` options installs requirements
needed for development, because CellRank is bundled with additional libraries.

Dependencies_
~~~~~~~~~~~~~
To efficiently compute the Schur decomposition for large cell numbers, we rely on `PETSc`_ and `SLEPc`_ which can
sometimes be a bit tricky to install. On Ubuntu 18.04, try the following::

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

During installation of petsc, petsc4py, slepc, and slepc4py the following
error might appear several times::

    ERROR: Failed building wheel for [insert package name here]

but this doesn't matter if the installer finally tells you::

    Successfully installed [insert package name here]

On Mac OS, install `MPICH`_ as a message passing interface and then proceed as above, using either pip or the
installation instructions given on the `PETSc`_ and `SLEPc`_ websites. The `SLEPc`_ homepage even offers a video tutorial
explaining the installation.

Note that this is only relevant for the :class:`cellrank.tl.GPCCA` estimator and only for large cell numbers.
For small cell numbers (<10k), you can use `method='brandts'` when computing the Schur decomposition.

Jupyter Notebook
~~~~~~~~~~~~~~~~

To run the tutorials in a notebook locally, please install::

   pip install notebook

and run ``jupyter notebook`` in the terminal. If you get the error ``Not a directory: 'xdg-settings'``,
use ``jupyter notebook --no-browser`` instead and open the url manually (or use this
`bugfix <https://github.com/jupyter/notebook/issues/3746#issuecomment-444957821>`_). Alternatively,
you can run all tutorials interactively directly in your browser using `binder`_. Just click the
binder button at the top of each tutorial.


If you run into issues, feel free to open a `GitHub issue`_ or send us an `email <mailto:info@cellrank.org>`_ .


.. _`Miniconda`: https://conda.pydata.org/miniconda.html
.. _`GitHub issue`: https://github.com/theislab/cellrank/issues/new
.. _`binder`: https://mybinder.org/
.. _`SLEPc`: https://slepc.upv.es/
.. _`PETSc`: https://www.mcs.anl.gov/petsc/
.. _`MPICH`: https://www.mpich.org/
