from typing import Any

from typing_extensions import Literal

from anndata import AnnData

from cellrank.ul._docs import d
from cellrank.tl._utils import cyto_trace
from cellrank.tl.kernels._pseudotime_kernel import PseudotimeKernel


@d.dedent
class CytoTRACEKernel(PseudotimeKernel):
    """
    Kernel which computes directed transition probabilities based on a KNN graph and the CytoTRACE score [Cyto20]_.

    The KNN graph contains information about the (undirected) connectivities among cells, reflecting their similarity.
    CytoTRACE can be used to estimate cellular plasticity and in turn, a pseudotemporal ordering of cells from more
    plactic to less plastic states. This kernel internally uses the `PseudotimeKernel` to direct the KNN graph on
    the basis of the CytoTRACE-derived pseudotime.

    %(density_correction)s

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :paramref:`adata` ``.obs`` where the pseudotime is stored.
    compute_cond_num
        Whether to compute condition number of the transition matrix. Note that this might be costly,
        since it does not use sparse implementation.
    kwargs
        Keyword arguments for :func:`cellrank.tl.cyto_trace`.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
        **kwargs: Any,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key="gcs_pseudotime",
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
            **kwargs,
        )
        self._time_key = "gcs_pseudotime"

    def _read_from_adata(
        self,
        time_key: str,
        layer: str = "Ms",
        aggregation: Literal["mean", "median", "hmean", "gmean"] = "mean",
        use_raw: bool = True,
        **kwargs: Any,
    ):
        cyto_trace(self.adata, layer=layer, aggregation=aggregation, use_raw=use_raw)

        super()._read_from_adata(time_key=time_key, **kwargs)
