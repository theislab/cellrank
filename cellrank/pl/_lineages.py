from typing import Any, Union, Optional, Sequence
from typing_extensions import Literal

from cellrank.pl._utils import AnnData
from cellrank.estimators import CFLARE
from cellrank._utils._key import Key
from cellrank._utils._docs import d

__all__ = ["lineages"]


@d.dedent
def lineages(
    adata: AnnData,
    lineages: Optional[Union[str, Sequence[str]]] = None,
    backward: bool = False,
    color: Optional[str] = None,
    mode: Literal["embedding", "time"] = "embedding",
    time_key: str = "latent_time",
    **kwargs: Any,
) -> None:
    """
    Plot lineages.

    For each lineage, we show all cells in an embedding (default is UMAP) and color them by their probability of
    belonging to this lineage. For cells that are already committed, this probability will be one for  their respective
    lineage and zero otherwise. For naive cells, these probabilities will be more balanced, reflecting
    the fact that naive cells have the potential to develop towards multiple endpoints.

    Parameters
    ----------
    %(adata)s
    lineages
        Plot only these lineages. If `None`, plot all lineages.
    %(backward)s
    color
        If given, plot cluster annotations left of the lineage probabilities.
    %(time_mode)s
    time_key
        Key in ``adata.obs`` where the pseudotime is stored.
    %(basis)s
    kwargs
        Keyword arguments for :meth:`cellrank.estimators.GPCCA.plot_absorption_probabilities`.

    Returns
    -------
    %(just_plots)s
    """

    # use CFLARE (no macrostates required)
    mc = CFLARE.from_adata(adata, obsp_key=Key.uns.kernel(backward))
    if mc.absorption_probabilities is None:
        raise RuntimeError("Compute absorption probabilities first.")

    # plot using the MC object
    color = kwargs.pop("cluster_key", color)
    mc.plot_absorption_probabilities(
        states=lineages,
        color=color,
        mode=mode,
        time_key=time_key,
        **kwargs,
    )
