from typing import Any, Optional

from abc import ABC, abstractmethod

from anndata import AnnData


class AnnDataMixin(ABC):
    """TODO."""

    def _read_from_adata(self, adata: AnnData, **kwargs: Any) -> bool:
        """TODO."""
        return True

    @property
    @abstractmethod
    def adata(self) -> AnnData:
        """TODO."""

    @adata.setter
    @abstractmethod
    def adata(self, adata: Optional[AnnData]) -> None:
        """TODO."""

    @abstractmethod
    def __len__(self) -> int:
        """TODO."""

    @abstractmethod
    def to_adata(self) -> AnnData:
        """TODO."""

    @classmethod
    @abstractmethod
    def from_adata(cls, adata: AnnData, **kwargs: Any) -> "AnnDataMixin":
        """TODO."""
        obj = cls(adata, **kwargs)
        obj._read_from_adata(adata)

        return obj
