from typing import Any, Union, Optional, Sequence
from typing_extensions import Literal, Protocol

import scvelo as scv
from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import _complex_warning

import numpy as np


class VectorPlotterProtocol(Protocol):  # noqa: D101
    @property
    def adata(self) -> AnnData:  # noqa: D102
        ...


class VectorPlotter:
    """Class that allow plotting of eigen- or Schur vectors in an embedding."""

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)

    @d.get_sections(base="plot_vectors", sections=["Parameters", "Returns"])
    @d.dedent
    def _plot_vectors(
        self: VectorPlotterProtocol,
        _which: Literal["eigen", "schur"],
        _vectors: np.ndarray,
        use: Optional[Union[int, Sequence[int]]] = None,
        abs_value: bool = False,
        cluster_key: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        """
        Plot vectors in an embedding.

        Parameters
        ----------
        use
            Which or how many vectors are to be plotted.
        abs_value
            Whether to take the absolute value before plotting.
        cluster_key
            Key in :attr:`anndata.AnnData.obs` for plotting categorical observations.
        %(basis)s
        kwargs
            Keyword arguments for :func:`scvelo.pl.scatter`.

        Returns
        -------
        %(just_plots)s
        """

        if _which not in ("eigen", "schur"):
            raise ValueError(
                f"Invalid kind `{_which!r}`. Valid options are: `eigen` or `schur`."
            )
        if _which == "schur":
            is_schur = True
            name = "Schur "
        else:
            is_schur = False
            name = "eigen"

        # check whether dimensions are consistent
        if self.adata.n_obs != _vectors.shape[0]:
            raise ValueError(
                f"Number of cells ({self.adata.n_obs}) is inconsistent with the first"
                f"dimension of vectors ({_vectors.shape[0]})."
            )

        if use is None:
            # fmt: off
            m = (
                getattr(self, "eigendecomposition").get("eigengap", _vectors.shape[1]) + 1
                if not is_schur and getattr(self, "eigendecomposition", None) is not None
                else _vectors.shape[1]
            )
            use = list(range(is_schur, m))
            # fmt: on
        elif isinstance(use, int):
            use = list(range(is_schur, use + is_schur))
        elif not isinstance(use, (np.ndarray, Sequence)):
            raise TypeError(
                f"Expected `use` to be `int` or a `Sequence`, found `{type(use).__name__!r}`."
            )
        use = list(use)
        if not use:
            raise ValueError("No vectors have been selected.")

        if max(use) >= _vectors.shape[1]:
            raise ValueError(
                f"Maximum specified {name}vector ({max(use)}) is larger "
                f"than the number of computed {name}vectors ({_vectors.shape[1]})."
            )
        V_ = _vectors[:, use]

        if is_schur:
            title = [f"{name}vector {i}" for i in use]
        else:
            D = kwargs.pop("D")
            V_ = _complex_warning(V_, use, use_imag=kwargs.pop("use_imag", False))
            title = [fr"$\lambda_{i}$={ev:.02f}" for i, ev in zip(use, D[use])]

        if abs_value:
            V_ = np.abs(V_)

        color = list(V_.T)
        if cluster_key is not None:
            color = [cluster_key] + color
        cmap = kwargs.pop("cmap", "viridis")

        logg.debug(f"Plotting `{use}` {name}vectors")

        scv.pl.scatter(self.adata, color=color, title=title, cmap=cmap, **kwargs)
