from typing import Any, Dict, Tuple, Union, Mapping, Optional

from tqdm.auto import tqdm
from typing_extensions import Literal

from anndata import AnnData

import numpy as np
import pandas as pd
from scipy.sparse import bmat, spdiags, spmatrix

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.tl._utils import _normalize

_error = None
try:
    import wot

    from cellrank.tl.kernels import ExperimentalTimeKernel as Kernel
except ImportError as e:
    from cellrank.external.kernels._import_error_kernel import ErroredKernel as Kernel

    _error = e
    wot = None


class nstr(str):  # used for params+cache (precomputed cost matrices)
    """String class that is not equal to any other string."""

    def __eq__(self, other: str) -> bool:
        return False


@d.dedent
class WOTKernel(Kernel, error=_error):
    """
    Waddington optimal transport kernel from [Schiebinger19]_.

    This class requires the `wot` package, which can be installed `pip install wot`.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    time_key
        Key in :attr:`adata` ``.obs`` where experimental time is stored.
        The experimental time can be of either of a numeric or an ordered categorical type.
    %(cond_num)s
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        time_key: str = "exp_time",
        compute_cond_num: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=time_key,
            compute_cond_num=compute_cond_num,
        )

        # WOT's requirements
        cats = self.experimental_time.cat.categories
        self._exp_time = self.experimental_time.cat.rename_categories(
            dict(zip(cats, map(float, cats)))
        )
        self.adata.obs[self._time_key] = self.experimental_time
        self._growth_rates = None

    def compute_transition_matrix(
        self,
        cost_matrices: Optional[
            Union[str, Mapping[Tuple[float, float], np.ndarray]]
        ] = None,
        solver: Literal["fixed_iters", "duality_gap"] = "duality_gap",
        growth_rate_key: Optional[str] = None,
        **kwargs: Any,
    ) -> "WOTKernel":
        """
        Compute transition matrix using Waddington OT [Schiebinger19]_.

        Parameters
        ----------
        cost_matrices
            Cost matrices for each consecutive time pair.
            If a :class:`str`, it specifies a key in :attr:`adata` ``.layers``: or :attr:`adata` ``.obsm``
            containing cell features that are used to compute cost matrices.
            If `None`, use `WOT`'s default.
        solver
            Which solver to use.
        growth_rate_key
            Key in :attr:`adata` ``.obs`` where initial cell growth rates are stored.
        kwargs
            Keyword arguments for OT configuration.

        Returns
        -------
        :class:`cellrank.external.kernels.WOTKernel`
            Makes :attr:`transition_matrix` and :attr:`growth_rates` available.
        """
        # disallow these params (because they e.g. subset the data)
        _ = kwargs.pop("cell_day_filter", None)
        _ = kwargs.pop("covariate_field", None)
        _ = kwargs.pop("ncounts", None)
        _ = kwargs.pop("ncells", None)
        kwargs.setdefault("epsilon", 0.05)
        kwargs.setdefault("lambda1", 1)
        kwargs.setdefault("lamdda2", 50)
        kwargs["growth_iters"] = max(kwargs.get("growth_iters"), 1)

        start = logg.info(
            "Computing transition matrix using Waddington Optimal Transport"
        )

        cost_matrices, cmat_param = self._generate_cost_matrices(cost_matrices)
        if self._reuse_cache(
            {
                "cost_matrices": cmat_param,
                "solver": solver,
                "growth_rate_key": growth_rate_key,
                **kwargs,
            },
            time=start,
        ):
            return self

        tmat = self._compute_pairwise_tmats(
            cost_matrices=cost_matrices,
            solver=solver,
            growth_rate_field=growth_rate_key,
            **kwargs,
        )
        tmat = self._restich_tmats(tmat)

        self._compute_transition_matrix(
            matrix=tmat,
            density_normalize=False,
            check_irreducibility=False,
        )

        logg.info("    Finish", time=start)

        return self

    def _compute_pairwise_tmats(
        self,
        cost_matrices: Optional[
            Union[str, Mapping[Tuple[float, float], np.ndarray]]
        ] = None,
        solver: Literal["fixed_iters", "duality_gap"] = "duality_gap",
        growth_rate_field: Optional[str] = None,
        **kwargs,
    ) -> Dict[Tuple[float, float], AnnData]:
        self._ot_model = wot.ot.OTModel(
            self.adata,
            day_field=self._time_key,
            covariate_field=None,
            growth_rate_field=growth_rate_field,
            solver=solver,
            **kwargs,
        )

        tmats: Dict[Tuple[float, float], AnnData] = {}
        start = logg.info(
            f"Computing transport maps for `{len(cost_matrices)}` time pairs"
        )
        for tpair, cost_matrix in tqdm(cost_matrices.items(), unit="time pair"):
            tmat: AnnData = self._ot_model.compute_transport_map(
                *tpair, cost_matrix=cost_matrix
            )
            nans = int(np.sum(np.isnan(tmat.X)))
            if nans:
                raise ValueError(
                    f"Encountered `{nans}` NaN values for time pair `{tpair}`."
                )
            tmat.X = _normalize(tmat.X)
            tmats[tpair] = tmat

        logg.info("    Finish", time=start)

        return tmats

    def _restich_tmats(self, tmaps: Mapping[Union[float, float], AnnData]) -> spmatrix:
        blocks = [[None] * (len(tmaps) + 1) for _ in range(len(tmaps) + 1)]
        nrows, ncols = 0, 0
        obs_names, growth_rates = [], []

        for i, tmap in enumerate(tmaps.values()):
            blocks[i][i + 1] = tmap.X
            nrows += tmap.n_obs
            ncols += tmap.n_vars
            obs_names.extend(tmap.obs_names)
            growth_rates.append(tmap.obs)
        obs_names.extend(tmap.var_names)

        n = self.adata.n_obs - nrows
        blocks[-1][-1] = spdiags([1] * n, 0, n, n)
        n = blocks[0][1].shape[0]
        blocks[0][0] = spdiags([], 0, n, n)

        tmp = AnnData(bmat(blocks, format="csr"))
        tmp.obs_names = obs_names
        tmp.var_names = obs_names
        tmp = tmp[self.adata.obs_names, :][:, self.adata.obs_names]

        self._growth_rates = pd.merge(
            tmp.obs,
            pd.concat(growth_rates),
            left_index=True,
            right_index=True,
            how="left",
        )

        return tmp.X

    def _generate_cost_matrices(
        self,
        cost_matrices: Optional[
            Union[str, Mapping[Tuple[float, float], np.ndarray]]
        ] = None,
    ) -> Tuple[Mapping[Tuple[float, float], Optional[np.ndarray]], str]:
        timepoints = self.experimental_time.cat.categories
        timepoints = list(zip(timepoints[:-1], timepoints[1:]))

        if cost_matrices is None:
            logg.debug("Using default cost matrices")
            return {tpair: None for tpair in timepoints}, "default"

        if isinstance(cost_matrices, dict):
            logg.debug("Using precomputed cost matrices")

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
            return cmats, nstr("precomputed")

        if isinstance(cost_matrices, str):
            logg.debug(f"Computing cost matrices using `{cost_matrices!r}` key")
            if cost_matrices == "X":
                cost_matrices = None

            try:
                features = self.adata._get_X(layer=cost_matrices)
                modifier = "layer"
            except KeyError:
                try:
                    features = self.adata.obsm[cost_matrices]
                    modifier = "obsm"
                except KeyError:
                    raise KeyError(
                        f"Unable to find key `{cost_matrices!r}` in `adata.layers` or `adata.obsm`."
                    ) from None

            cmats = {}
            for tpair in timepoints:
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
