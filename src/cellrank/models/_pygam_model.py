import collections
import copy
import enum
import types
import warnings
from typing import Any, Literal, Mapping, Optional, Union

from pygam import GAM as pGAM
from pygam import (
    ExpectileGAM,
    GammaGAM,
    InvGaussGAM,
    LinearGAM,
    LogisticGAM,
    PoissonGAM,
    s,
)

import numpy as np

from anndata import AnnData

from cellrank import logging as logg
from cellrank._utils._docs import d
from cellrank._utils._enum import ModeEnum
from cellrank._utils._utils import _filter_kwargs
from cellrank.models import BaseModel

__all__ = ["GAM"]


class GamLinkFunction(ModeEnum):
    IDENTITY = enum.auto()
    LOGIT = enum.auto()
    INVERSE = enum.auto()
    LOG = enum.auto()
    INV_SQUARED = enum.auto()


class GamDistribution(ModeEnum):
    NORMAL = enum.auto()
    BINOMIAL = enum.auto()
    POISSON = enum.auto()
    GAMMA = enum.auto()
    GAUSSIAN = enum.auto()
    INV_GAUSS = enum.auto()


_gams = collections.defaultdict(
    lambda: pGAM,
    {
        (GamDistribution.NORMAL, GamLinkFunction.IDENTITY): LinearGAM,
        (GamDistribution.BINOMIAL, GamLinkFunction.LOGIT): LogisticGAM,
        (GamDistribution.POISSON, GamLinkFunction.LOG): PoissonGAM,
        (GamDistribution.GAMMA, GamLinkFunction.LOG): GammaGAM,
        (GamDistribution.INV_GAUSS, GamLinkFunction.LOG): InvGaussGAM,
    },
)


@d.dedent
class GAM(BaseModel):
    """Fit Generalized Additive Models (GAMs) using :mod:`pygam`.

    Parameters
    ----------
    %(adata)s
    n_knots
        Number of knots.
    spline_order
        Order of the splines, e.g., :math:`3` for cubic splines.
    distribution
        Name of the distribution. Available distributions can be found
        `here <https://pygam.readthedocs.io/en/latest/notebooks/tour_of_pygam.html#Distribution:>`__.
    link
        Name of the link function. Available link functions can be found
        `here <https://pygam.readthedocs.io/en/latest/notebooks/tour_of_pygam.html#Link-function:>`__.
    max_iter
        Maximum number of iterations for optimization.
    expectile
        Expectile for :class:`~pygam.pygam.ExpectileGAM`. This forces the distribution to be ``'normal'``
        and link function to ``'identity'``. Must be in :math:`(0, 1)`.
    grid
        Whether to perform a grid search. Keys correspond to a parameter names and values to range to be searched.
        If ``'default'``, use the default grid. If :obj:`None`, don't perform a grid search.
    spline_kwargs
        Keyword arguments for :func:`~pygam.terms.s`.
    kwargs
        Keyword arguments for the :class:`~pygam.pygam.GAM`.
    """

    def __init__(
        self,
        adata: AnnData,
        n_knots: Optional[int] = 6,
        spline_order: int = 3,
        distribution: Literal["normal", "binomial", "poisson", "gamma", "gaussian", "inv_gauss"] = "gamma",
        link: Literal["identity", "logit", "inverse", "log", "inv_squared"] = "log",
        max_iter: int = 2000,
        expectile: Optional[float] = None,
        grid: Optional[Union[str, Mapping[str, Any]]] = None,
        spline_kwargs: Mapping[str, Any] = types.MappingProxyType({}),
        **kwargs: Any,
    ):
        term = s(
            0,
            spline_order=spline_order,
            n_splines=n_knots,
            **_filter_kwargs(s, **{**{"lam": 3, "penalties": ["derivative", "l2"]}, **spline_kwargs}),
        )
        link = GamLinkFunction(link)
        distribution = GamDistribution(distribution)
        if distribution == GamDistribution.GAUSSIAN:
            distribution = GamDistribution.NORMAL

        if expectile is not None:
            if not (0 < expectile < 1):
                raise ValueError(f"Expected `expectile` to be in `(0, 1)`, found `{expectile}`.")
            if distribution != "normal" or link != "identity":
                logg.warning(
                    f"Expectile GAM works only with `normal` distribution and `identity` link function,"
                    f"found `{distribution!r}` distribution and {link!r} link functions."
                )
            model = ExpectileGAM(term, expectile=expectile, max_iter=max_iter, verbose=False, **kwargs)
        else:
            # doing it like this ensure that user can specify scale
            gam = _gams[distribution, link]

            filtered_kwargs = _filter_kwargs(gam.__init__, **kwargs)
            if len(kwargs) != len(filtered_kwargs):
                raise TypeError(
                    f"Invalid arguments `{list(set(kwargs) - set(filtered_kwargs))}` " f"for `{type(gam).__name__!r}`."
                )

            filtered_kwargs["link"] = link
            filtered_kwargs["distribution"] = distribution

            model = gam(
                term,
                max_iter=max_iter,
                verbose=False,
                **_filter_kwargs(gam.__init__, **filtered_kwargs),
            )
        super().__init__(adata, model=model)

        if grid is None:
            self._grid = None
        elif isinstance(grid, dict):
            self._grid = copy.copy(grid)
        elif isinstance(grid, str):
            self._grid = object() if grid == "default" else None
        else:
            raise TypeError(f"Expected `grid` to be `dict`, `str` or `None`, found `{type(grid).__name__!r}`.")

    @d.dedent
    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> "GAM":
        """%(base_model_fit.full_desc)s

        Parameters
        ----------
        %(base_model_fit.parameters)s

        Returns
        -------
        Fits the model and returns self.
        """  # noqa

        super().fit(x, y, w, **kwargs)

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=DeprecationWarning,
                message=".* is a deprecated alias for the builtin",
            )
            if self._grid is not None:
                # use default search
                grid = {} if not isinstance(self._grid, dict) else self._grid
                try:
                    self.model.gridsearch(
                        self.x,
                        self.y,
                        weights=self.w,
                        keep_best=True,
                        progress=False,
                        **grid,
                        **kwargs,
                    )
                    return self
                except Exception as e:  # noqa: BLE001
                    # workaround for: https://github.com/dswah/pyGAM/issues/273
                    self.model.fit(self.x, self.y, weights=self.w, **kwargs)
                    logg.error(f"Grid search failed, reason: `{e}`. Fitting with default values")

            try:
                self.model.fit(self.x, self.y, weights=self.w, **kwargs)
                return self
            except Exception as e:  # noqa: BLE001
                raise RuntimeError(
                    f"Unable to fit `{type(self).__name__}` for gene "
                    f"`{self._gene!r}` in lineage `{self._lineage!r}`. Reason: `{e}`"
                ) from e

    @d.dedent
    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs: Any,
    ) -> np.ndarray:
        """%(base_model_predict.full_desc)s

        Parameters
        ----------
        %(base_model_predict.parameters)s

        Returns
        -------
        %(base_model_predict.returns)s
        """  # noqa

        x_test = self._check(key_added, x_test)

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=DeprecationWarning,
                message=".* is a deprecated alias for the builtin",
            )
            self._y_test = self.model.predict(x_test, **kwargs)
        self._y_test = np.squeeze(self._y_test).astype(self._dtype)

        return self.y_test

    @d.dedent
    def confidence_interval(self, x_test: Optional[np.ndarray] = None, **kwargs: Any) -> np.ndarray:
        """%(base_model_ci.summary)s

        Parameters
        ----------
        %(base_model_ci.parameters)s

        Returns
        -------
        %(base_model_ci.returns)s
        """  # noqa

        x_test = self._check("_x_test", x_test)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=DeprecationWarning,
                message=".* is a deprecated alias for the builtin",
            )
            self._conf_int = self.model.confidence_intervals(x_test, **kwargs).astype(self._dtype)

        return self.conf_int

    @d.dedent
    def copy(self) -> "BaseModel":
        """%(copy)s"""  # noqa
        res = GAM(self.adata)
        self._shallowcopy_attributes(res)
        res._grid = copy.deepcopy(self._grid)
        return res
