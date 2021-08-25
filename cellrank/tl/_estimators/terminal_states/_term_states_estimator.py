from typing import Any, Union, Optional

from abc import ABC, abstractmethod

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
        self._term_states_probs: Optional[pd.Series] = None
        self._term_states_colors: Optional[np.ndarray] = None
        # TODO: set_terminal_states + rename_terminal_states

    def to_adata(self) -> None:
        super().to_adata()

        key = Key.obs.term_states(self.backward)
        # TODO: set_or_debug
        if self.terminal_states is not None:
            self.adata.obs[key] = self.terminal_states
        if self.terminal_states_probabilities is not None:
            self.adata.obs[Key.obs.probs(key)] = self.terminal_states_probabilities
        if self._term_states_colors is not None:
            self.adata.uns[Key.uns.colors(key)] = self.terminal_states

    def fit(self, *args: Any, **kwargs: Any) -> None:
        # TODO: implement me
        return NotImplemented

    @abstractmethod
    def compute_terminal_states(self, *args: Any, **kwargs: Any) -> None:
        pass

    @property
    def terminal_states(self) -> Optional[pd.Series]:
        """TODO."""
        return self._term_states

    @property
    def terminal_states_probabilities(self) -> Optional[pd.Series]:
        """TODO."""
        return self._term_states_probs
