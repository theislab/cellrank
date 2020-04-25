# -*- coding: utf-8 -*-
from typing import Optional

from cellrank.tools._estimators._base_estimator import BaseEstimator


class GPCCA(BaseEstimator):
    def _read_from_adata(
        self, g2m_key: Optional[str] = None, s_key: Optional[str] = None, **kwargs
    ) -> None:
        raise NotImplementedError()

    def copy(self) -> "GPCCA":
        raise NotImplementedError()
