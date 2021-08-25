from typing import Union, Optional

from abc import ABC

from anndata import AnnData
from cellrank.tl._estimators import BaseEstimator
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl._estimators.mixins._constants import Key

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix


class TermStatesEstimator(BaseEstimator, ABC):
    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        key: Optional[str] = None,
        obsp_key: Optional[str] = None,
    ):
        super().__init__(obj=obj, key=key, obsp_key=obsp_key)
        self._term_states: Optional[pd.Series] = None
        self._term_states_colors: Optional[np.ndarray] = None
        # TODO: abstract compute_term_states
        # TODO: implement fit here
        # TODO: set_terminal_states + rename_terminal_states

    def to_adata(self) -> None:
        super().to_adata()

        key = Key.obs.term_states(self.backward)
        if self.terminal_states is not None:
            self.adata.obs[key] = self.terminal_states
        if self._term_states_colors is not None:
            self.adata.uns[Key.uns.colors(key)] = self.terminal_states

    @property
    def terminal_states(self) -> Optional[pd.Series]:
        """TODO."""
        return self._term_states
