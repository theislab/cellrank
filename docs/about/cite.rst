Citing CellRank
===============
If you find CellRank useful for your research, please consider citing our work as follows: If you are using
CellRank's _VelocityKernel_ with classical RNA velocity, cite :cite:`lange:22` as:

.. code-block:: bibtex

    @article{lange:22,
        title = {CellRank for directed single-cell fate mapping},
        author = {Lange, Marius and Bergen, Volker and Klein, Michal and Setty, Manu and
                  Reuter, Bernhard and Bakhti, Mostafa and Lickert, Heiko and
                  Ansari, Meshal and Schniering, Janine and Schiller, Herbert B. and
                  Pe'er, Dana and Theis, Fabian J.},
        journal = {Nat. Methods},
        year = {2022},
        doi = {10.1038/s41592-021-01346-6},
        publisher = {Nature Publishing Group}
    }

If you are using the _PseudotimeKernel_, _CytoTRACEKernel_, `RealTimeKernel`, or _VelocityKernel_ with RNA velocity inferred
from metabolic labeling data using our approach, cite :cite:`weiler:23` as:

.. code-block:: bibtex

    @article{weiler:23,
        title = {Unified fate mapping in multiview single-cell data},
        author = {Weiler, Philipp and Lange, Marius and Klein, Michal and Pe{\textquotesingle}er, Dana and Theis, Fabian},
        doi = {10.1101/2023.07.19.549685},
        url = {https://doi.org/10.1101/2023.07.19.549685},
        year = {2023},
        publisher = {Cold Spring Harbor Laboratory},
    }

In addition, if you use the :class:`~cellrank.estimators.GPCCA` estimator to compute initial, terminal or intermediate
states, you are using the `pyGPCCA package <https://github.com/msmdev/pyGPCCA>`_ :cite:`reuter:22` under the hood,
which implements the Generalized Perron Cluster Cluster Analysis (GPCCA) algorithm. Thus, additionally to CellRank,
please cite GPCCA :cite:`reuter:19` as:

.. code-block:: bibtex

    @article{reuter:19,
        author = {Reuter,Bernhard  and Fackeldey,Konstantin  and Weber,Marcus },
        title = {Generalized Markov modeling of nonreversible molecular kinetics},
        journal = {The Journal of Chemical Physics},
        volume = {150},
        number = {17},
        pages = {174103},
        year = {2019},
        doi = {10.1063/1.5064530},
    }
