Installation
============
CellRank requires Python3 version >= 3.6 to run. We recommend using Miniconda_.

GitHub
~~~~~~
Install CellRank by calling::

    pip install git+https://github.com/theislab/cellrank

Development Version
~~~~~~~~~~~~~~~~~~~
To stay up-to-date with the newest version, run::

    git clone https://github.com/theislab/cellrank && cd cellrank
    pip install -e .

``-e`` stands for ``--editable`` and makes sure that your environment is updated
when you pull new changes from github.

Dependencies
^^^^^^^^^^^^
To efficiently compute the schur decomposition for large cell numbers, we rely on `SLEPSc`_ wich can
sometimes be a bit tricky to install. On Ubuntu 18.04, try the following::

    # update
    sudo apt-get update
    sudo apt-get upgrade

    # install a message passing interface and mpi4py
    sudo apt-get install libopenmpi-dev
    pip install --user mpi4py

    # install petsc and and petsc4py
    pip install --user petsc
    pip install --user petsc4py

    # install slepsc and slepsc4py
    pip install --user slepc slepc4py

During installation of petsc, petsc4py, selpc, and slepc4py the following
error might appear several times::

ERROR: Failed building wheel for [insert package name here]

but this doesn't matter if the installer finally tells you::

Successfully installed [insert package name here]

Jupyter Notebook
^^^^^^^^^^^^^^^^

To run the tutorials in a notebook locally, please install::

   pip install notebook

and run ``jupyter notebook`` in the terminal. If you get the error ``Not a directory: 'xdg-settings'``,
use ``jupyter notebook --no-browser`` instead and open the url manually (or use this
`bugfix <https://github.com/jupyter/notebook/issues/3746#issuecomment-444957821>`_). Alternatively,
you can run all tutorials interactively directly in your browser using `binder`_. Just click the
binder button at the top of each tutorial.


If you run into issues, feel free to open a `GitHub issue`_ or send us an `email <mailto:info@cellrank.org>`_ .


.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _`Github issue`: https://github.com/theislab/cellrank/issues/new
.. _`binder`: https://mybinder.org/
.. _`SLEPSc`: https://slepc.upv.es/
