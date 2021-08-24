from typing import Any, Union, Optional, Sequence
from typing_extensions import Literal, Protocol

import scvelo as scv
from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import _complex_warning

import numpy as np


class VectorPlottableProtocol(Protocol):
    @property
    def adata(self) -> AnnData:
        ...


class VectorPlottable:
    @d.get_sections(base="plot_vectors", sections=["Parameters", "Returns"])
    @d.dedent
    def _plot_vectors(
        self: VectorPlottableProtocol,
        which: Literal["eigen", "schur"],
        vectors: np.ndarray,
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

        if which not in ("eigen", "schur"):
            raise ValueError(
                f"Invalid kind `{which!r}`. Valid options are `eigen` or `schur`."
            )
        if which == "schur":
            is_schur = True
            name = "Schur "
        else:
            is_schur = False
            name = "eigen"

        # check whether dimensions are consistent
        if self.adata.n_obs != vectors.shape[0]:
            raise ValueError(
                f"Number of cells ({self.adata.n_obs}) is inconsistent with the first"
                f"dimension of vectors ({vectors.shape[0]})."
            )

        if use is None:
            # fmt: off
            m = (
                getattr(self, "eigendecomposition").get("eigengap", vectors.shape[1]) + 1
                if not is_schur and getattr(self, "eigendecomposition", None) is not None
                else vectors.shape[1]
            )
            use = list(range(is_schur, m))
            # fmt: on
        elif isinstance(use, int):
            use = list(range(is_schur, use + is_schur))
        elif not isinstance(use, (range, np.ndarray, Sequence)):
            raise TypeError(
                f"Argument `use` must be either `int`, `range` or a `sequence`, "
                f"found `{type(use).__name__}`."
            )
        else:
            if not all(map(lambda u: isinstance(u, int), use)):
                raise TypeError("Not all values in `use` argument are integers.")
        use = list(use)
        if not use:
            raise ValueError("Nothing to plot.")

        muse = max(use)
        if muse >= vectors.shape[1]:
            raise ValueError(
                f"Maximum specified {name}vector ({muse}) is larger "
                f"than the number of computed {name}vectors ({vectors.shape[1]})."
            )
        V_ = vectors[:, use]

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
