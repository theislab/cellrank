from typing import Any, Dict, Type, Tuple, Union, Literal, Mapping, Optional, Sequence

from types import MappingProxyType
from moscot.base.output import BaseSolverOutput
from moscot.problems.base import CompoundProblem

import scanpy as sc
import anndata
from cellrank.kernels._transport_map_kernel import TransportMapKernel


class MoscotKernel(TransportMapKernel):
    """Kernel for all moscot problems with a transport matrix as output."""

    def __init__(
        self, *args: Any, problem: Optional[Type[CompoundProblem]] = None, **kwargs: Any
    ):
        super().__init__(*args, **kwargs)
        self._problem = problem

    @classmethod
    def load(
        cls,
        problem: Type[CompoundProblem],
        backward: bool = False,
        **kwargs: Any,
    ) -> "MoscotKernel":
        """Load from a moscot model."""
        return MoscotKernel(
            adata=problem.adata,
            time_key=problem.temporal_key,
            backward=backward,
            problem=problem,
            **kwargs,
        )

    def compute_transition_matrix(
        self,
        tmaps: Optional[Dict[Tuple[float, float], Type[BaseSolverOutput]]] = None,
        threshold: Optional[str] = None,
        self_transitions: Union[
            Literal["uniform", "diagonal", "connectivities", "all"],
            Sequence[float],
        ] = "connectivities",
        conn_weight: Optional[float] = None,
        conn_kwargs: Mapping[str, Any] = MappingProxyType({}),
    ) -> "MoscotKernel":
        """Write transport matrices to kernel."""
        if tmaps is None:
            if self.problem is None:
                raise ValueError(
                    "If `MoscotKernel.problem` is `None`, `tmaps` must be provided."
                )
            self._moscot_to_cr_tmaps({p.solution for p in self.problem.problems})
        else:
            self._tmaps = self._moscot_to_cr_tmaps(tmaps)

        if threshold is not None:
            self._threshold_transport_maps(
                self.transport_maps, threshold=threshold, copy=False
            )

        tmap = self._restich_tmaps(
            self.transport_maps,
            self_transitions=self_transitions,
            conn_weight=conn_weight,
            **conn_kwargs,
        )

        self.transition_matrix = tmap.X
        return self

    def _compute_tmap(self):
        raise NotImplementedError()

    def _moscot_to_cr_tmaps(
        self, outputs: Dict[Tuple[float, float], Type[BaseSolverOutput]]
    ) -> Dict[Tuple[float, float], anndata.AnnData]:
        tmats = {}
        for (t1, t2), output in outputs.items():
            src_mask = self.adata.obs[self._time_key] == t1
            tgt_mask = self.adata.obs[self._time_key] == t2
            tmats[t1, t2] = anndata.AnnData(
                output.transport_matrix,
                obs=self.adata[src_mask].obs,
                var=self.adata[tgt_mask].obs,
            )

    @property
    def problem(self) -> Type[CompoundProblem]:
        """Return the moscot problem."""
        return self._problem
