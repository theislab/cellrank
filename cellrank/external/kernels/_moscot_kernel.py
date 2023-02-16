from typing import Any, Type, Tuple, Union, Literal, Mapping, Optional, Sequence

from types import MappingProxyType

import scanpy as sc
from anndata import AnnData
from cellrank import logging as logg
from cellrank._utils._docs import d
from cellrank._utils._utils import _maybe_subset_hvgs
from cellrank.external.kernels._utils import MarkerGenes
from cellrank.kernels.utils._tmat_flow import Numeric_t
from cellrank.kernels._transport_map_kernel import Pair_t, Threshold_t, SelfTransitions

import numpy as np
import pandas as pd

__all__ = ["MoscotKernel"]

_error = None
try:
    import moscot
    from moscot.problems.base import CompoundProblem

    from cellrank.kernels import BaseTransportMapKernel as Kernel
except ImportError as e:
    from cellrank.external.kernels._import_error_kernel import ErroredKernel as Kernel

    _error = e
    moscot = None


# TODO(michalk8): refactor me properly using TransportMapKernel + update snippet
@d.dedent
class MoscotKernel(Kernel, error=_error):
    """Kernel to interface with moscot problems."""

    @classmethod
    def load(
        cls,
        problem: Type[CompoundProblem],
        time_points: Optional[Tuple[float, float]] = None,
        threshold: Optional[float] = None,
        self_transitions: Union[
            Literal["uniform", "diagonal", "connectivities", "all"],
            Sequence[Numeric_t],
        ] = SelfTransitions.CONNECTIVITIES,
        conn_weight: Optional[float] = None,
        conn_kwargs: Mapping[str, Any] = MappingProxyType({}),
        **kwargs: Any,
    ) -> "MoscotKernel":
        """Load moscot solution into CellRank kernel."""
        kernel = MoscotKernel(
            adata=problem.adata, time_key=problem.temporal_key, backward=False, **kwargs
        )

        if time_points is None:
            time_points = list(problem)
        if threshold is None:
            kernel.transport_maps = [
                problem[key].solution.transport_matrix for key in time_points
            ]
        else:
            kernel.transport_maps = [
                problem[key].solution.sparsify(threshold=threshold)
                for key in time_points
            ]

        return kernel
