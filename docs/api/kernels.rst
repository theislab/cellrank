.. module:: cellrank.kernels

Kernels
=======
Kernels compute cell-cell transition probabilities based on various input data modalities, including molecular
similarity, RNA velocity :cite:`manno:18`, experimental time points, and many more. They come with methods
for high-level, qualitative visualization, including vector field- or random walk plots. For quantitative analysis of
kernel-computed transition matrices, we recommend using :doc:`estimators <estimators>`.

We don't use the term "kernel" in the way it is used in mathematics, but rather colloquially, to refer to a class that
takes in multi-view single-cell data and outputs a cell-cell transition matrix.

.. currentmodule:: cellrank

.. autosummary::
    :toctree: _autosummary/kernels

    kernels.VelocityKernel
    kernels.ConnectivityKernel
    kernels.PseudotimeKernel
    kernels.CytoTRACEKernel
    kernels.TransportMapKernel
    kernels.PrecomputedKernel
