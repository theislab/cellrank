from typing import Any

from abc import ABC, abstractmethod

from anndata import AnnData


class AnnDataMixin(ABC):
    """TODO."""

    def __init__(self, adata: AnnData, **kwargs: Any):
        super().__init__(**kwargs)

        self._adata = adata

    @property
    def adata(self) -> AnnData:
        """TODO."""
        return self._adata

    def __len__(self) -> int:
        return len(self.adata)

    @abstractmethod
    def to_adata(self) -> None:
        """TODO."""

    @classmethod
    @abstractmethod
    def from_adata(cls, adata: AnnData) -> "AnnDataMixin":
        """TODO."""
