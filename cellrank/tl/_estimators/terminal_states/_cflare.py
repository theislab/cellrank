from typing import Any

from cellrank.tl._estimators.mixins import EigenMixin
from cellrank.tl._estimators.terminal_states._term_states_estimator import (
    TermStatesEstimator,
)


class CFLARE(TermStatesEstimator, EigenMixin):
    def fit(self, *args: Any, **kwargs: Any) -> None:
        return NotImplemented
