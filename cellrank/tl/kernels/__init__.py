# -*- coding: utf-8 -*-
from cellrank.tl.kernels._base_kernel import Kernel, Constant
from cellrank.tl.kernels._palantir_kernel import PalantirKernel
from cellrank.tl.kernels._velocity_kernel import VelocityKernel
from cellrank.tl.kernels._velocity_schemes import (
    CosineScheme,
    DotProductScheme,
    CorrelationScheme,
    SimilaritySchemeABC,
)
from cellrank.tl.kernels._precomputed_kernel import PrecomputedKernel
from cellrank.tl.kernels._connectivity_kernel import ConnectivityKernel
