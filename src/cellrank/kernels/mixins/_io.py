import contextlib
import pathlib
import pickle
from typing import Optional, Protocol, Tuple, Union

from anndata import AnnData

from cellrank import logging as logg

__all__ = ["IOMixin"]


class IOMixinProtocol(Protocol):
    @property
    def shape(self) -> Tuple[int, ...]:
        ...

    @property
    def adata(self) -> AnnData:
        ...

    @adata.setter
    def adata(self, adata: AnnData) -> None:
        ...


class IOMixin:
    """Mixin that allows for serialization from/to files using :mod:`pickle`."""

    @property
    @contextlib.contextmanager
    def _remove_adata(self) -> None:
        """Temporarily remove :attr:`adata`, if present."""
        adata = getattr(self, "adata", None)

        try:
            if adata is not None:
                self.adata = None
            yield
        finally:
            if adata is not None:
                self.adata = adata

    def write(
        self,
        fname: Union[str, pathlib.Path],
        write_adata: bool = True,
    ) -> None:
        """Serialize self to a file using :mod:`pickle`.

        Parameters
        ----------
        fname
            Path where to save the object.
        write_adata
            Whether to save :attr:`adata` object.

        Returns
        -------
        Nothing, just writes itself to a file.
        """
        logg.info(f"Writing `{self}` to `{fname}`")

        if write_adata:
            with open(fname, "wb") as fout:
                pickle.dump(self, fout)
        else:
            with self._remove_adata, open(fname, "wb") as fout:
                pickle.dump(self, fout)

    @staticmethod
    def read(fname: Union[str, pathlib.Path], adata: Optional[AnnData] = None, copy: bool = False) -> "IOMixin":
        """De-serialize self from a file.

        Parameters
        ----------
        fname
            Path from which to read the object.
        adata
            :class:`~anndata.AnnData` object to assign to the saved object.
            Only used when the saved object has :attr:`adata` and it was saved without it.
        copy
            Whether to copy ``adata`` before assigning it. If ``adata`` is a view, it is always copied.

        Returns
        -------
        The de-serialized object.
        """
        with open(fname, "rb") as fin:
            obj: IOMixinProtocol = pickle.load(fin)

        if hasattr(obj, "adata"):
            if isinstance(obj.adata, AnnData):
                if adata is not None:
                    logg.warning("Ignoring supplied `adata` object because it is already present")
                return obj

            if not isinstance(adata, AnnData):
                raise TypeError(
                    "This object was saved without its `adata` object. " "Please supply one as `adata=...`."
                )

            if obj.shape[0] != len(adata):
                raise ValueError(f"Expected `adata` to be of length `{len(adata)}`, found `{obj.shape[0]}`.")
            if copy or adata.is_view:
                adata = adata.copy()

            obj.adata = adata
            return obj

        return obj
