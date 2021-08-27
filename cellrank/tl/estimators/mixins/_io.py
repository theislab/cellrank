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

    @property
    @contextmanager
    def _remove_adata(self) -> None:
        adata = getattr(self, "_adata", None)

        try:
            if adata is not None:
                self._adata = None
            yield
        finally:
            if adata is not None:
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

        if keep_adata:
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
            TODO.
        copy
            TODO.

        Returns
        -------
        The deserialized object.
        """

        with open(fname, "rb") as fin:
            obj = pickle.load(fin)

        if hasattr(obj, "_adata"):
            if not isinstance(adata, AnnData):
                raise TypeError(
                    "This object was saved without its `anndata.AnnData` object. "
                    "Please supply one as `adata=...`."
                )
            try:
                assert len(obj) == len(adata), "TODO."
                # if hasattr(obj, "obs_names"):
                #    adata = adata[obj.obs_names]
                if copy or adata.is_view:
                    adata = adata.copy()
            except AssertionError as e:
                raise ValueError(e) from None
            except TypeError as e:
                raise AttributeError(e) from None
            except KeyError as e:
                raise KeyError(e) from None

            obj._adata = adata
            return obj

        return obj
