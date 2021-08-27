from typing import Any, Optional
from typing_extensions import Protocol

from anndata import AnnData
from cellrank import logging as logg

import pandas as pd
from scipy.stats import ranksums
from pandas.api.types import infer_dtype, is_categorical_dtype


class CCDetectorProtocol(Protocol):
    @property
    def adata(self) -> AnnData:
        ...


class CCDetectorMixin:
    def __init__(
        self: CCDetectorProtocol,
        g2m_key: Optional[str] = None,
        s_key: Optional[str] = None,
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        # fmt: off
        try:
            self._g2m_score: Optional[pd.Series] = None if g2m_key is None else self.adata.obs[g2m_key].copy()
            self._s_score: Optional[pd.Series] = None if s_key is None else self.adata.obs[g2m_key].copy()
        except KeyError as e:
            raise KeyError(f"Unable to find cell-cycle score in `adata.obs[{str(e)}]`.") from None
        # fmt: on

    def _detect_cc_stages(self, categories: pd.Series, p_thresh: float = 1e-15) -> None:
        """
        Detect cell-cycle driven start or endpoints.

        Parameters
        ----------
        categories
            Categorical values for which to determine whether they are cell-cycle-driver or not.
        p_thresh
            p-value threshold for the rank-sum test for the group to be considered cell-cycle driven.

        Returns
        -------
        Nothing, just warns if a group in ``categories`` is cell-cycle driven.
        """
        if not is_categorical_dtype(categories):
            raise TypeError(
                f"Object must be `categorical`, found `{infer_dtype(categories).__name__}`."
            )

        groups = list(categories.cat.categories)
        if not groups:
            return

        scores = []
        if self._g2m_score is not None:
            scores.append(self._g2m_score)
        if self._s_score is not None:
            scores.append(self._s_score)
        if not scores:
            return

        for group in groups:
            mask = categories == group
            for score in scores:
                a, b = score[mask], score[~mask]
                statistic, pvalue = ranksums(a, b)
                if statistic > 0 and pvalue < p_thresh:
                    logg.warning(f"Group `{group!r}` appears to be cell-cycle driven")
                    break
