# -*- coding: utf-8 -*-
"""Module containing :mod:`pygam` model implementation."""
from copy import copy as _copy
from copy import deepcopy
from types import MappingProxyType
from typing import Mapping, Optional
from collections import defaultdict

import numpy as np
from pygam import GAM as pGAM
from pygam import (
    GammaGAM,
    LinearGAM,
    PoissonGAM,
    InvGaussGAM,
    LogisticGAM,
    ExpectileGAM,
    s,
)

from cellrank import logging as logg
from cellrank.ul._docs import d
from cellrank.ul.models import BaseModel
from cellrank.tl._constants import ModeEnum
from cellrank.tl.kernels._utils import _filter_kwargs
from cellrank.ul.models._base_model import AnnData

_r_lib = None
_r_lib_name = None


class GamLinkFunction(ModeEnum):  # noqa
    IDENTITY = "identity"
    LOGIT = "logit"
    INV = "inverse"
    LOG = "log"
    INV_SQUARED = "inverse-squared"


class GamDistribution(ModeEnum):  # noqa
    NORMAL = "normal"
    BINOMIAL = "binomial"
    POISSON = "poisson"
    GAMMA = "gamma"
    GAUSS = "gaussian"
    INV_GAUSS = "inv_gauss"


_gams = defaultdict(
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
    """
    Fit Generalized Additive Models (GAMs) using package :mod:`pygam`.

    Parameters
    ----------
    %(adata)s
    n_splines
        Number of splines.
    spline_order
        Order of the splines.
    distribution
        Name of the distribution. Available distributions can be found
        `here <https://pygam.readthedocs.io/en/latest/notebooks/tour_of_pygam.html#Distribution:>`_.
    link
        Name of the link function. Available link functions can be found
        `here <https://pygam.readthedocs.io/en/latest/notebooks/tour_of_pygam.html#Link-function:>`_.
    max_iter
        Maximum number of iterations for optimization.
    expectile
        Expectile for :class:`pygam.pygam.ExpectileGAM`. This forces the distribution to be `'normal'`
        and link function to `'identity'`. Must be in `(0, 1)`.
    use_default_conf_int
        Whether to use :meth:`default_confidence_interval` to calculate the confidence interval or
        use the :paramref:`model`'s method.
    grid
        Whether to perform a grid search. Keys correspond to a parameter names and values to range to be searched.
        If an empty :class:`dict`, don't perform a grid search. If `None`, use a default grid.
    spline_kwargs
        Keyword arguments for :class:`pygam.s`.
    **kwargs
        Keyword arguments for :class:`pygam.pygam.GAM`.
    """

    def __init__(
        self,
        adata: AnnData,
        n_splines: Optional[int] = 10,
        spline_order: int = 3,
        distribution: str = "gamma",
        link: str = "log",
        max_iter: int = 2000,
        expectile: Optional[float] = None,
        use_default_conf_int: bool = False,
        grid: Optional[Mapping] = None,
        spline_kwargs: Mapping = MappingProxyType({}),
        **kwargs,
    ):
        term = s(
            0,
            spline_order=spline_order,
            n_splines=n_splines,
            penalties=["derivative", "l2"],
            **_filter_kwargs(s, **{**{"lam": 3}, **spline_kwargs}),
        )
        link = GamLinkFunction(link)
        distribution = GamDistribution(distribution)
        if distribution == GamDistribution.GAUSS:
            distribution = GamDistribution.NORMAL

        if expectile is not None:
            if not (0 < expectile < 1):
                raise ValueError(
                    f"Expected `expectile` to be in `(0, 1)`, found `{expectile}`."
                )
            if distribution != "normal" or link != "identity":
                logg.warning(
                    f"Expectile GAM works only with `normal` distribution and `identity` link function,"
                    f"found `{distribution!r}` distribution and {link!r} link functions."
                )
            model = ExpectileGAM(
                term, expectile=expectile, max_iter=max_iter, verbose=False, **kwargs
            )
        else:
            gam = _gams[
                distribution, link
            ]  # doing it like this ensure that user can specify scale
            kwargs["link"] = link.s
            kwargs["distribution"] = distribution.s
            model = gam(
                term,
                max_iter=max_iter,
                verbose=False,
                **_filter_kwargs(gam.__init__, **kwargs),
            )
        super().__init__(adata, model=model)
        self._use_default_conf_int = use_default_conf_int

        if grid is None:
            self._grid = None
        elif isinstance(grid, dict):
            self._grid = _copy(grid)
        elif isinstance(grid, str):
            self._grid = object() if grid == "default" else None
        else:
            raise TypeError(
                f"Expected `grid` to be `dict`, `str` or `None`, found `{type(grid).__name__!r}`."
            )

    @d.dedent
    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "GAM":
        """
        %(base_model_fit.full_desc)s

        Parameters
        ----------
        %(base_model_fit.parameters)s

        Returns
        -------
        :class:`cellrank.ul.models.GAM`
            Fits the model and returns self.
        """  # noqa

        super().fit(x, y, w, **kwargs)

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
            except Exception as e:
                # workaround for: https://github.com/dswah/pyGAM/issues/273
                self.model.fit(self.x, self.y, weights=self.w, **kwargs)
                logg.error(
                    f"Grid search failed, reason: `{e}`. Fitting with default values"
                )

        try:
            self.model.fit(self.x, self.y, weights=self.w, **kwargs)
            return self
        except Exception as e:
            raise RuntimeError(
                f"Unable to fit `{type(self).__name__}` for gene "
                f"`{self._gene!r}` in lineage `{self._lineage!r}`. Reason: `{e}`"
            ) from e

    @d.dedent
    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: Optional[str] = "_x_test",
        **kwargs,
    ) -> np.ndarray:
        """
        %(base_model_predict.full_desc)s

        Parameters
        ----------
        %(base_model_predict.parameters)s

        Returns
        -------
        %(base_model_predict.returns)s
        """  # noqa

        x_test = self._check(key_added, x_test)

        self._y_test = self.model.predict(x_test, **kwargs)
        self._y_test = np.squeeze(self._y_test)

        return self.y_test

    @d.dedent
    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs
    ) -> np.ndarray:
        """
        %(base_model_ci.summary)s

        Parameters
        ----------
        %(base_model_ci.parameters)s

        Returns
        -------
        %(base_model_ci.returns)s
        """  # noqa

        x_test = self._check("_x_test", x_test)
        if self._use_default_conf_int:
            self._conf_int = self.default_confidence_interval(x_test=x_test, **kwargs)
        else:
            self._conf_int = self.model.confidence_intervals(x_test, **kwargs)

        return self.conf_int

    @d.dedent
    def copy(self) -> "BaseModel":
        """%(copy)s"""  # noqa
        res = GAM(self.adata)

        res._use_default_conf_int = self._use_default_conf_int
        res._grid = deepcopy(self._grid)
        res._model = deepcopy(self.model)

        return res
