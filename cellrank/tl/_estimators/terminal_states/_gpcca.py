from typing import Any, Union, Optional

from anndata import AnnData
from cellrank.tl._estimators.mixins import EigenMixin, SchurMixin, LinDriversMixin
from cellrank.tl.kernels._base_kernel import KernelExpression
from cellrank.tl._estimators.terminal_states._term_states_estimator import (
    TermStatesEstimator,
)

import numpy as np
from scipy.sparse import spmatrix


class GPCCA(TermStatesEstimator, LinDriversMixin, SchurMixin, EigenMixin):
    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        key: Optional[str] = None,
        obsp_key: Optional[str] = None,
    ):
        super().__init__(obj=obj, key=key, obsp_key=obsp_key)

    def fit(self, *args: Any, **kwargs: Any) -> None:
        return NotImplemented
