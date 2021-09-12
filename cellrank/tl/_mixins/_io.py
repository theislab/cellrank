from typing import Any, Tuple, Union, Optional
from typing_extensions import Protocol

import pickle
from pathlib import Path
from contextlib import contextmanager

from anndata import AnnData
from cellrank import logging as logg


class IOMixinProtocol(Protocol):  # noqa: D101
    @property
    def shape(self) -> Tuple[int, ...]:  # noqa: D102
        ...

    @property
    def adata(self) -> AnnData:  # noqa: D102
        ...

    @adata.setter
    def adata(self, adata: AnnData) -> None:
        ...


class IOMixin:
    """Mixin that allows for serialization from/to files using :mod:`pickle`."""

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)

    @property
    @contextmanager
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
        fname: Union[str, Path],
        write_adata: bool = True,
        ext: Optional[str] = "pickle",
    ) -> None:
        """
        Serialize self to a file.

        Parameters
        ----------
        fname
            Filename where to save the object.
        write_adata
            Whether to save :attr:`adata` object or not, if present.
        ext
            Filename extension to use. If `None`, don't append any extension.

        Returns
        -------
        Nothing, just writes itself to a file using :mod:`pickle`.
        """

        fname = str(fname)
        if ext is not None:
            if not ext.startswith("."):
                ext = "." + ext
            if not fname.endswith(ext):
                fname += ext

        logg.info(f"Writing `{self}` to `{fname}`")

        if write_adata:
            with open(fname, "wb") as fout:
                pickle.dump(self, fout)
            return

        with self._remove_adata:
            with open(fname, "wb") as fout:
                pickle.dump(self, fout)

    @staticmethod
    def read(
        fname: Union[str, Path], adata: Optional[AnnData] = None, copy: bool = False
    ) -> "IOMixin":
        """
        Deserialize self from a file.

        Parameters
        ----------
        fname
            Filename from which to read the object.
        adata
            :class:`anndata.AnnData` object to assign to the saved object.
            Only used when the saved object has :attr:`adata` and it was saved without it.
        copy
            Whether to copy ``adata`` before assigning it or not. If ``adata`` is a view, it is always copied.

        Returns
        -------
        The deserialized object.
        """

        with open(fname, "rb") as fin:
            obj: IOMixinProtocol = pickle.load(fin)

        if hasattr(obj, "adata"):
            if isinstance(obj.adata, AnnData):
                if adata is not None:
                    logg.warning(
                        "Ignoring supplied `adata` object because it is already present"
                    )
                return obj

            if not isinstance(adata, AnnData):
                raise TypeError(
                    "This object was saved without its `adata` object. "
                    "Please supply one as `adata=...`."
                )

            if obj.shape[0] != len(adata):
                raise ValueError(
                    f"Expected `adata` to be of length `{len(adata)}`, found `{obj.shape[0]}`."
                )
            if copy or adata.is_view:
                adata = adata.copy()

            obj.adata = adata
            return obj

        return obj
