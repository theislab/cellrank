from typing import Any

from abc import ABC, abstractmethod

from anndata import AnnData


class AnnDataMixin(ABC):
    """TODO."""

    def __init__(self, adata: AnnData, **kwargs: Any):
        super().__init__(**kwargs)

        self._adata = adata
        self._n_obs = len(adata)

    @property
    def adata(self) -> AnnData:
        """TODO."""
        return self._adata

    def __len__(self) -> int:
        return self._n_obs

    @abstractmethod
    def to_adata(self) -> AnnData:
        """TODO."""

    def _read_from_adata(self, adata: AnnData, **kwargs: Any) -> bool:
        """TODO."""
        return True

    @classmethod
    @abstractmethod
    def from_adata(cls, adata: AnnData, **kwargs: Any) -> "AnnDataMixin":
        """TODO."""
