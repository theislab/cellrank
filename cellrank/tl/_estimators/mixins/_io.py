from typing import Any, Union, Optional

import pickle
from pathlib import Path
from contextlib import contextmanager

from anndata import AnnData
from cellrank import logging as logg


class IOMixin:
    """TODO."""

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)

    @contextmanager
    def _maybe_remove_adata(self, remove: bool) -> None:
        adata = getattr(self, "_adata", None)
        should_remove = remove and adata is not None
        setter = False

        try:
            if should_remove:
                try:
                    # invoke setter if possible
                    self.adata = None
                    setter = True
                except AttributeError:
                    self._adata = None
            yield
        finally:
            if should_remove:
                if setter:
                    self.adata = adata
                else:
                    self._adata = adata

    def write(
        self,
        fname: Union[str, Path],
        keep_adata: bool = True,
        ext: Optional[str] = "pickle",
    ) -> None:
        """
        Serialize self to a file.

        Parameters
        ----------
        fname
            Filename where to save the object.
        keep_adata
            TODO.
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

        with self._maybe_remove_adata(not keep_adata):
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
            TODO.
        copy
            TODO.

        Returns
        -------
        The deserialized object.
        """

        with open(fname, "rb") as fin:
            obj = pickle.load(fin)

        if adata is None:
            return obj
        if not isinstance(adata, AnnData):
            raise TypeError("TODO")

        if hasattr(obj, "_adata"):
            try:
                assert len(obj) == len(adata)
                if hasattr(obj, "obs_names"):
                    adata = adata[obj.obs_names]
                if copy or adata.is_view:
                    adata = adata.copy()
            except AssertionError:
                raise ValueError("TODO") from None
            except TypeError:
                raise AttributeError("TODO") from None
            except KeyError:
                raise KeyError("TODO") from None
            obj._adata = adata
            return obj

        logg.warning("TODO")
        return obj
