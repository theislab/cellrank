from typing import Any, Union, Optional

from anndata import AnnData
from cellrank.tl import Lineage
from cellrank.tl._estimators.mixins import EigenMixin, SchurMixin, LinDriversMixin
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl._estimators.terminal_states._term_states_estimator import (
    TermStatesEstimator,
)

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix


class GPCCA(TermStatesEstimator, LinDriversMixin, SchurMixin, EigenMixin):
    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        obsp_key: Optional[str] = None,
        **kwargs: Any,
    ):
        super().__init__(obj=obj, obsp_key=obsp_key, **kwargs)

        self._coarse_init_dist: Optional[pd.Series] = None
        self._coarse_stat_dist: Optional[pd.Series] = None
        self._coarse_tmat: Optional[pd.DataFrame] = None

        self._macrostates: Optional[Lineage] = None
        # TODO: do we need this?
        self._macrostates_probs: Optional[pd.Series] = None
        self._macrostates_colors: Optional[np.ndarray] = None

    def compute_terminal_states(self, *args: Any, **kwargs: Any) -> None:
        pass

    def fit(self, *args: Any, **kwargs: Any) -> None:
        # TOOO: call super + optionally abs prob?
        return NotImplemented

    @property
    def macrostates(self) -> Optional[Lineage]:
        """TODO."""
        return self._macrostates

    @property
    def coarse_initial_distribution(self) -> Optional[pd.Series]:
        """TODO."""
        return self._coarse_init_dist

    @property
    def coarse_stationary_distribution(self) -> Optional[pd.Series]:
        """TODO."""
        return self._coarse_stat_dist

    @property
    def coarse_T(self) -> Optional[pd.DataFrame]:
        """TODO."""
        return self._coarse_tmat