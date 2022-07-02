.. _kernels:

Kernels
=======
TODO(Marius1311): be more descriptive/verbose.
Kernels are used to estimate cell-to-cell transitions.

.. module:: cellrank.kernels
.. currentmodule:: cellrank

.. autosummary::
    :toctree: _autosummary/kernels

    kernels.VelocityKernel
    kernels.ConnectivityKernel
    kernels.PseudotimeKernel
    kernels.CytoTRACEKernel
    kernels.PrecomputedKernel

External Kernels
----------------
.. module:: cellrank.external
.. currentmodule:: cellrank

.. autosummary::
    :toctree: _autosummary/kernels

    external.kernels.StationaryOTKernel
    external.kernels.WOTKernel
