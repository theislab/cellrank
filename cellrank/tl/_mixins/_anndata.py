from typing import Any, Optional

from abc import ABC, abstractmethod

from anndata import AnnData
from cellrank.ul._docs import d


class AnnDataMixin(ABC):
    """Mixin that allows for serialization from/to :class:`anndata.AnnData`."""

    @abstractmethod
    @d.dedent
    def _read_from_adata(self, adata: AnnData, **kwargs: Any) -> bool:
        """
        Populate attributes of self from :class:`anndata.AnnData`.

        Parameters
        ----------
        %(adata)s
        kwargs
            Additional keyword arguments.

        Returns
        -------
        `True` if the deserialization was successful, otherwise `False`.
        """
        return True

    @property
    @abstractmethod
    def adata(self) -> AnnData:
        """Annotated data object."""

    @adata.setter
    @abstractmethod
    def adata(self, adata: Optional[AnnData]) -> None:
        ...

    @abstractmethod
    def __len__(self) -> int:
        ...

    @abstractmethod
    @d.get_full_description(base="to_adata")
    def to_adata(self, **kwargs: Any) -> AnnData:
        """Serialize self to :class:`anndata.Anndata`."""

    @classmethod
    @d.get_full_description(base="from_adata")
    @d.get_sections(base="from_adata", sections=["Returns"])
    @d.dedent
    def from_adata(cls, adata: AnnData, **kwargs: Any) -> "AnnDataMixin":
        """
        Deserialize self from :class:`anndata.AnnData`.

        Parameters
        ----------
        %(adata)s
        kwargs
            Additional keyword arguments.

        Returns
        -------
        The deserialized object.
        """
        obj = cls(adata, **kwargs)
        obj._read_from_adata(adata)

        return obj
