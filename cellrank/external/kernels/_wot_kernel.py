from typing import Any, Dict, Tuple, Union, Mapping, Optional, Sequence
from typing_extensions import Literal

from types import MappingProxyType
from tqdm.auto import tqdm

import scanpy as sc
from anndata import AnnData
from cellrank import logging as logg
from cellrank._utils._docs import d
from cellrank._utils._utils import _maybe_subset_hvgs
from cellrank.external.kernels._utils import MarkerGenes
from cellrank.kernels.utils._tmat_flow import Numeric_t
from cellrank.kernels._transport_map_kernel import Pair_t, Threshold_t, SelfTransitions

import numpy as np
import pandas as pd

__all__ = ["WOTKernel"]

_error = None
try:
    import wot

    from cellrank.kernels import TransportMapKernel as Kernel
except ImportError as e:
    from cellrank.external.kernels._import_error_kernel import ErroredKernel as Kernel

    _error = e
    wot = None


# TODO(michalk8): refactor me properly using TransportMapKernel
@d.dedent
class WOTKernel(Kernel, error=_error):
    """
    Waddington optimal transport kernel from :cite:`schiebinger:19`.

    This class requires the `wot` package, which can be installed as
    `pip install git+https://github.com/broadinstitute/wot`.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :attr:`anndata.AnnData.obs` where experimental time is stored.
        The experimental time can be of either of a numeric or an ordered categorical type.
    kwargs
        Keyword arguments for the parent class.

    Examples
    --------
    Workflow::

        # import packages, load data
        import scanpy as sc
        import cellrank as cr
        adata = cr.datasets.lung()

        # filter, normalize and annotate highly variable genes
        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)

        # estimate proliferation and apoptosis from gene sets (see. e.g. WOT tutorial for example lists)
        proliferation_genes = ...
        apoptosis_genes = ...
        sc.tl.score_genes(adata, gene_list=proliferation_genes, score_name='proliferation')
        sc.tl.score_genes(adata, gene_list=apoptosis_genes, score_name='apoptosis')

        # initialize kernel, estimate initial growth rate based on scores from above
        from cellrank.external.kernels import WOTKernel
        ot = WOTKernel(adata, time_key='day')
        ot.compute_initial_growth_rates(proliferation_key='proliferation',
                                        apoptosis_key='apoptosis',
                                        key_added='initial_growth_rates')

        # compute transport maps, aggregate into one single transition matrix
        ot.compute_transition_matrix(growth_rate_key='initial_growth_rates', growth_iters=3)
    """

    __import_error_message__ = (
        "Unable to import the kernel. Please install `wot` first as "
        "`pip install git+https://github.com/broadinstitute/wot`."
    )

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        time_key: str = "exp_time",
        **kwargs: Any,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=time_key,
            **kwargs,
        )

        # WOT's requirements
        cats = self.experimental_time.cat.categories
        self._exp_time = self.experimental_time.cat.rename_categories(
            dict(zip(cats, map(float, cats)))
        )
        self.adata.obs[self._time_key] = self.experimental_time
        self._growth_rates = None

    def compute_initial_growth_rates(
        self,
        proliferation_key: Optional[str] = None,
        apoptosis_key: Optional[str] = None,
        organism: Optional[Literal["human", "mouse"]] = None,
        beta_min: float = 0.3,
        beta_max: float = 1.7,
        delta_min: float = 0.3,
        delta_max: float = 1.7,
        key_added: Optional[str] = None,
        **kwargs: Any,
    ) -> Optional[pd.Series]:
        r"""
        Estimate initial growth rates using a birth-death process as described in :cite:`schiebinger:19`.

        The doubling time is defined as :math:`\frac{\ln 2}{\beta - \delta}` (similarly defined for half-time).
        The logistic function is used to transform the birth/death rate scores and to smoothly interpolate between
        specified minimum and maximum birth/death rates.

        Parameters
        ----------
        proliferation_key
            Key in :attr:`anndata.AnnData.obs` where the birth rate score is saved.
        apoptosis_key
            Key in :attr:`anndata.AnnData.obs` where the death rate score is saved.
        organism
            Organism for which to calculate the birth/death scores, if they cannot be found in :attr:`adata`.
            In this case, :func:`scanpy.tl.score_genes` is used to calculate the scores based on an organism-dependent
            set of marker genes for proliferation and apoptosis.
        beta_min
            Minimum birth rate.
        beta_max
            Maximum birth rate.
        delta_min
            Minimum death rate.
        delta_max
            Maximum death rate.
        key_added
            Key in :attr:`adata` ``.obs`` where to add the estimated growth rates. If `None`, just return them.
        kwargs
            Keyword arguments for :func:`scanpy.tl.score_genes`. Only used when ``proliferation_key``
            or ``apoptosis_key`` cannot be found in :attr:`adata.AnnData.obs`.

        Returns
        -------
        :class:`pandas.Series`
            The estimated initial growth rates if ``key_added = None``, otherwise `None`.

        Notes
        -----
        If you don't have access to proliferation/apoptosis gene sets, you can use the ones defined in :mod:`cellrank`
        for a specific organism. Alternatively, you can also use WOT without an estimate of initial growth rates. In
        that case, make sure to use several iterations in
        :meth:`cellrank.external.kernels.WOTKernel.compute_transition_matrix` by increasing the ``growth_iters``
        parameter. A value around 3 works well in most cases.

        The markers used here were taken from the following sources:

            - human, proliferation - :cite:`tirosh:16:science`.
            - human, apoptosis - `Hallmark Apoptosis, MSigDB <https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_APOPTOSIS>`_.
            - mouse, proliferation - :cite:`tirosh:16:nature`.
            - mouse, apoptosis - `Hallmark P53 Pathway, MSigDB <https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_P53_PATHWAY>`_.

        For more information about WOT, see the official `tutorial <https://broadinstitute.github.io/wot/tutorial/>`_.
        """  # noqa: E501

        def get_scores(
            key: str,
            *,
            kind: Literal["proliferation", "apoptosis"],
            organism: Optional[Literal["human", "mouse"]],
        ) -> np.ndarray:
            try:
                return np.asarray(self.adata.obs[key])
            except KeyError:
                if organism is None:
                    raise KeyError(
                        f"Unable to find `{kind}` scores in `adata.obs[{kind!r}]`. Consider specifying `organism=...`."
                    ) from None

                logg.info(f"Computing `{kind}` scores")
                score_name = f"{kind}_score" if key is None else key

                sc.tl.score_genes(
                    self.adata,
                    gene_list=getattr(MarkerGenes, f"{kind}_markers")(organism),
                    score_name=score_name,
                    **kwargs,
                )

                return get_scores(score_name, kind=kind, organism=None)

        def logistic(x: np.ndarray, L: float, k: float, x0: float = 0) -> np.ndarray:
            return L / (1 + np.exp(-k * (x - x0)))

        def gen_logistic(
            p: np.ndarray, beta_max: float, beta_min: float, center: float, width: float
        ) -> np.ndarray:
            return beta_min + logistic(p, L=beta_max - beta_min, k=4 / width, x0=center)

        def beta(
            p: np.ndarray,
            beta_max: float = 1.7,
            beta_min: float = 0.3,
            center: float = 0.25,
        ) -> np.ndarray:
            return gen_logistic(p, beta_max, beta_min, center, width=0.5)

        def delta(
            a: np.ndarray,
            delta_max: float = 1.7,
            delta_min: float = 0.3,
            center: float = 0.1,
        ) -> np.ndarray:
            return gen_logistic(a, delta_max, delta_min, center, width=0.2)

        birth = beta(
            get_scores(proliferation_key, kind="proliferation", organism=organism),
            beta_min=beta_min,
            beta_max=beta_max,
        )
        death = delta(
            get_scores(apoptosis_key, kind="apoptosis", organism=organism),
            delta_min=delta_min,
            delta_max=delta_max,
        )
        gr = np.exp(birth - death)

        if key_added is None:
            return pd.Series(gr, index=self.adata.obs_names)
        self.adata.obs[key_added] = gr

    @d.dedent
    def compute_transition_matrix(
        self,
        cost_matrices: Optional[Union[str, Mapping[Pair_t, np.ndarray]]] = None,
        lambda1: float = 1,
        lambda2: float = 50,
        epsilon: float = 0.05,
        growth_iters: int = 1,
        solver: Literal["fixed_iters", "duality_gap"] = "duality_gap",
        growth_rate_key: Optional[str] = None,
        use_highly_variable: Optional[Union[str, bool]] = True,
        threshold: Optional[Threshold_t] = "auto",
        self_transitions: Union[
            Literal["uniform", "diagonal", "connectivities", "all"],
            Sequence[Union[Numeric_t]],
        ] = SelfTransitions.CONNECTIVITIES,
        conn_weight: Optional[float] = None,
        conn_kwargs: Mapping[str, Any] = MappingProxyType({}),
        **kwargs: Any,
    ) -> "WOTKernel":
        """
        Compute transition matrix using Waddington OT :cite:`schiebinger:19`.

        Computes transport maps linking together pairs of time points for time-series single cell data using unbalanced
        optimal transport, taking into account cell birth and death rates. From the sequence of transition maps linking
        pairs of sequential time points, we construct one large transition matrix which contains the normalized
        transport maps as blocks on the 1st upper diagonal.

        Parameters
        ----------
        cost_matrices
            Cost matrices for each consecutive pair of time points.
            If a :class:`str`, it specifies a key in :attr:`anndata.AnnData.layers` or :attr:`anndata.AnnData.obsm`.
            containing cell features that are used to compute cost matrices. If `None`, use `WOT`'s default, i.e.
            compute distances in PCA space derived from :attr:`anndata.AnnData.X` for each time point pair separately.
        lambda1
            Regularization parameter for the marginal constraint on :math:`p`, the transport map row sums.
            Smaller value is useful when precise information about the growth rate is not present.
        lambda2
            Regularization parameter for the marginal constraint on :math:`q`, the transport map column sums.
        epsilon
            Entropy regularization parameter. Larger value gives more entropic descendant distributions.
        growth_iters
            Number of iterations for growth rate estimates. If growth rates are not known, consider using more
            iterations.
        solver
            Which solver to use.
        growth_rate_key
            Key in :attr:`anndata.AnnData.obs` where initial cell growth rates are stored.
            See :meth:`compute_initial_growth_rates` on how to estimate them.
        use_highly_variable
            Key in :attr:`anndata.AnnData.var` where highly variable genes are stored.
            If `True`, use `'highly_variable'`. If `None`, use all genes.
        %(tmk_tmat.parameters)s
        kwargs
            Additional keyword arguments for optimal transport configuration.

        Returns
        -------
        Self and updates the following attributes:

            - :attr:`transition_matrix` - transition matrix.
            - :attr:`transport_maps` - transport maps between consecutive time points.
            - :attr:`growth_rates` - estimated growth rates.

        Also modifies :attr:`anndata.AnnData.obs` with the following key:

            - `'estimated_growth_rates'` - the final estimated growth rates.

        Notes
        -----
        For more information about WOT, see the official `tutorial <https://broadinstitute.github.io/wot/tutorial/>`_.
        """
        # disallow these params (because they e.g. subset the data)
        _ = kwargs.pop("cell_day_filter", None)
        _ = kwargs.pop("covariate_field", None)
        _ = kwargs.pop("ncounts", None)
        _ = kwargs.pop("ncells", None)
        _ = kwargs.pop("parameters", None)  # parameters file
        kwargs["lambda1"] = lambda1
        kwargs["lambda2"] = lambda2
        kwargs["epsilon"] = epsilon
        kwargs["growth_iters"] = max(growth_iters, 1)
        self_transitions = SelfTransitions(self_transitions)

        start = logg.info(
            "Computing transition matrix using Waddington optimal transport"
        )

        adata = _maybe_subset_hvgs(self.adata, use_highly_variable=use_highly_variable)
        cost_matrices, cmat_param = self._generate_cost_matrices(adata, cost_matrices)
        if self._reuse_cache(
            {
                "cost_matrices": cmat_param,
                "solver": solver,
                "growth_rate_key": growth_rate_key,
                "use_highly_variable": use_highly_variable,
                "self_transitions": self_transitions,
                "threshold": threshold,
                **kwargs,
            },
            time=start,
        ):
            return self

        logg.info(f"Computing transport maps for `{len(cost_matrices)}` time pairs")
        self._tmaps = self._compute_tmaps(
            adata,
            cost_matrices=cost_matrices,
            solver=solver,
            growth_rate_field=growth_rate_key,
            **kwargs,
        )
        if threshold is not None:
            self._threshold_transport_maps(self._tmaps, threshold=threshold, copy=False)

        tmap = self._restich_tmaps(
            self._tmaps,
            self_transitions=self_transitions,
            conn_weight=conn_weight,
            conn_kwargs=conn_kwargs,
        )
        self._growth_rates = tmap.obs

        self.transition_matrix = tmap.X
        self.adata.obs["estimated_growth_rates"] = self.growth_rates[f"g{growth_iters}"]

        logg.info("    Finish", time=start)

        return self

    def _compute_tmap(self, t1: Numeric_t, t2: Numeric_t, **kwargs: Any) -> AnnData:
        ot_model: wot.ot.OTModel = kwargs.pop("ot_model")
        cost_matrix: Optional[Union[str, np.ndarray]] = kwargs.pop("cost_matrix")

        tmap: Optional[AnnData] = ot_model.compute_transport_map(
            t1,
            t2,
            cost_matrix=cost_matrix,
        )
        if tmap is None:
            raise TypeError(
                f"Unable to compute transport map for time pair `{(t1, t2)}`. "
                f"Please ensure that `adata.obs[{self._time_key!r}]` has the correct dtype (float)."
            )
        tmap.X = tmap.X.astype(np.float64)
        nans = int(np.sum(~np.isfinite(tmap.X)))
        if nans:
            raise ValueError(
                f"Encountered `{nans}` non-finite values for time pair `{(t1, t2)}`."
            )

        return tmap

    def _compute_tmaps(
        self,
        adata: AnnData,
        cost_matrices: Optional[Mapping[Pair_t, Optional[np.ndarray]]] = None,
        solver: Literal["fixed_iters", "duality_gap"] = "duality_gap",
        growth_rate_field: Optional[str] = None,
        **kwargs: Any,
    ) -> Dict[Tuple[float, float], AnnData]:
        _ = wot.ot.OTModel(
            adata,
            day_field=self._time_key,
            covariate_field=None,
            growth_rate_field=growth_rate_field,
        )
        for k in kwargs:
            if k not in _.ot_config:
                raise TypeError(f"WOT got an unexpected keyword argument {k!r}.")

        ot_model = wot.ot.OTModel(
            adata,
            day_field=self._time_key,
            covariate_field=None,
            growth_rate_field=growth_rate_field,
            solver=solver,
            **kwargs,
        )

        tmaps: Dict[Tuple[float, float], AnnData] = {}
        for (t1, t2), cost_matrix in tqdm(cost_matrices.items(), unit="time pair"):
            tmaps[t1, t2] = self._compute_tmap(
                t1,
                t2,
                ot_model=ot_model,
                cost_matrix=cost_matrix,
            )

        return tmaps

    def _generate_cost_matrices(
        self,
        adata: AnnData,
        cost_matrices: Optional[Union[str, Mapping[Pair_t, np.ndarray]]] = None,
    ) -> Tuple[Mapping[Pair_t, Optional[np.ndarray]], str]:
        timepoints = self.experimental_time.cat.categories
        timepoints = list(zip(timepoints[:-1], timepoints[1:]))

        if cost_matrices is None:
            logg.info("Using default cost matrices")
            return {tpair: None for tpair in timepoints}, "default"

        if isinstance(cost_matrices, dict):
            logg.info("Using precomputed cost matrices")

            cmats = {}
            for tpair in timepoints:
                if tpair not in cost_matrices:
                    logg.warning(
                        f"Unable to find cost matrix for pair `{tpair}`. Using default"
                    )
                cmats[tpair] = cmat = cost_matrices.get(tpair, None)
                if cmat is not None:
                    n_start = len(np.where(self.experimental_time == tpair[0])[0])
                    n_end = len(np.where(self.experimental_time == tpair[1])[0])
                    try:
                        if cmat.shape != (n_start, n_end):
                            raise ValueError(
                                f"Expected cost matrix for time pair `{tpair}` to be "
                                f"of shape `{(n_start, n_end)}`, found `{cmat.shape}`."
                            )
                    except AttributeError:
                        logg.warning(
                            f"Unable to verify whether supplied cost matrix for time pair `{tpair}` "
                            f"has the correct shape `{(n_start, n_end)}`"
                        )

            # prevent equality comparison when comparing with cache
            return cmats, "precomputed"

        if isinstance(cost_matrices, str):
            logg.info(f"Computing cost matrices using `{cost_matrices!r}` key")
            if cost_matrices == "X":
                cost_matrices = None

            try:
                features = adata._get_X(layer=cost_matrices)
                modifier = "layer"
            except KeyError:
                try:
                    features = adata.obsm[cost_matrices]
                    modifier = "obsm"
                except KeyError:
                    raise KeyError(
                        f"Unable to find key `{cost_matrices!r}` in `adata.layers` or `adata.obsm`."
                    ) from None

            cmats = {}
            for tpair in tqdm(timepoints, unit="cost matrix"):
                start_ixs = np.where(self.experimental_time == tpair[0])[0]
                end_ixs = np.where(self.experimental_time == tpair[1])[0]

                # being sparse is handled in WOT's function below
                cmats[tpair] = wot.ot.OTModel.compute_default_cost_matrix(
                    features[start_ixs], features[end_ixs]
                )

            return cmats, f"{modifier}:{cost_matrices}"

        raise NotImplementedError(
            f"Specifying cost matrices as "
            f"`{type(cost_matrices).__name__}` is not yet implemented."
        )

    @property
    def growth_rates(self) -> Optional[pd.DataFrame]:
        """Estimated cell growth rates for each growth rate iteration."""
        return self._growth_rates

    def __invert__(self) -> "WOTKernel":
        wk = super().__invert__()
        # needed for WOT
        wk.adata.obs[wk._time_key] = wk.experimental_time
        return wk
