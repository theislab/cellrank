# -*- coding: utf-8 -*-
"""Module containing all models somehow interfacing R."""
from copy import copy
from typing import Any, Tuple, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm

from cellrank.ul._docs import d
from cellrank.ul.models import BaseModel
from cellrank.ul.models._utils import _find_knots, _get_offset
from cellrank.ul.models._base_model import AnnData

_r_lib = None
_r_lib_name = None


@d.dedent
class GAMR(BaseModel):
    """
    Wrapper around R's `mgcv <https://cran.r-project.org/web/packages/mgcv/>`_ or
    `gam <https://cran.r-project.org/web/packages/gam/>`_ package for fitting Generalized Additive Models (GAMs).

    Parameters
    ----------
    %(adata)s
    n_splines
        Number of splines for the GAM.
    distribution
        Distribution family in `rpy2.robjects.r`, such as `'gaussian'`, `'poisson'` or `'nb'`.
        If `'nb'`, we always use the data in :paramref:`adata` ``.raw``.
    **kwargs
        TODO.
    """  # noqa

    def __init__(
        self,
        adata: AnnData,
        n_splines: int = 5,
        distribution: str = "gaussian",
        offset: Optional[np.ndarray] = None,
        perform_import_check: bool = True,
        **kwargs,
    ):
        super().__init__(adata, model=None)
        self._n_splines = n_splines
        self._family = distribution
        self._formula = f"y ~ s(x0, bs='cr', k={self._n_splines})"
        self._design_mat = None
        self._original_offset = None  # for copy
        self._gam_kwargs = copy(kwargs)

        self._lib = None
        self._lib_name = None

        if perform_import_check:
            # it's a bit costly to import, copying just passes the reference
            self._lib, self._lib_name = _maybe_import_r_lib("mgcv")

        if distribution == "nb":
            if offset is None:
                offset = _get_offset(adata, use_raw=True)
            offset = np.asarray(offset, dtype=self._dtype)

            if offset.shape != (adata.n_obs,):
                raise ValueError(
                    f"Expected offset to be of shape `{(adata.n_obs,)}`, found `{offset.shape}`."
                )

            # TODO: it's still WIP
            self._formula += " + offset(offset)"
            self._original_offset = offset

    def prepare(
        self,
        *args,
        **kwargs,
    ) -> "GAMR":
        """TODO"""  # noqa

        if self._family == "nb":
            kwargs["use_raw"] = True
        _ = super().prepare(*args, **kwargs)

        use_ixs = self.w > 0
        self._x = self.x[use_ixs]
        self._y = self.y[use_ixs]
        self._w = self.w[use_ixs]

        self._design_mat = pd.DataFrame(
            np.c_[self.x, self.y],
            columns=["x0", "y"],
        )

        if self._original_offset is not None:
            mask = np.isin(self.adata.obs_names, self._obs_names[use_ixs])
            self._design_mat["offset"] = self._original_offset[mask]

        return self

    @d.dedent
    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "GAMR":
        """
        %(base_model_fit.full_desc)s

        Parameters
        ----------
        %(base_model_fit.parameters)s

        Returns
        -------
        :class:`cellrank.ul.models.GAMR`
            Fits the model and returns self. Updates the following fields by filtering out `0` weights :paramref:`w`:

                - :paramref:`x` - %(base_model_x.summary)s
                - :paramref:`y` - %(base_model_y.summary)s
                - :paramref:`w` - %(base_model_w.summary)s
        """  # noqa

        from rpy2.robjects import Formula, r, pandas2ri

        super().fit(x, y, w, **kwargs)

        pandas2ri.activate()

        family = getattr(r, self._family)

        self._model = self._lib.gam(
            Formula(self._formula),
            data=self._design_mat,
            family=family,
            weights=pd.Series(self.w),
            knots=pd.DataFrame(
                _find_knots(self.x, self._n_splines)
            ),  # needs to be a DataFrame
            **self._gam_kwargs,
        )

        pandas2ri.deactivate()

        return self

    def _get_x_test(
        self, x_test: Optional[np.ndarray] = None, key_added: str = "_x_test"
    ) -> pd.DataFrame:
        newdata = pd.DataFrame(self._check(key_added, x_test), columns=["x0"])
        if "offset" in self._design_mat:
            newdata["offset"] = np.mean(self._design_mat["offset"])

        return newdata

    @d.dedent
    def predict(
        self, x_test: Optional[np.ndarray] = None, key_added: str = "_x_test", **kwargs
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

        from rpy2 import robjects
        from rpy2.robjects import pandas2ri

        if self.model is None:
            raise RuntimeError(
                "Trying to call an uninitialized model. To initialize it, run `.fit()` first."
            )
        if self._lib is None:
            raise RuntimeError(
                f"Unable to run prediction, R package `{self._lib_name!r}` is not imported."
            )

        newdata = self._get_x_test(x_test, key_added=key_added)

        pandas2ri.activate()
        self._y_test = (
            np.array(
                robjects.r.predict(
                    self.model,
                    newdata=pandas2ri.py2rpy(newdata),
                )
            )
            .squeeze()
            .astype(self._dtype)
        )

        pandas2ri.deactivate()

        return self.y_test

    @d.dedent
    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, level: float = 0.95, **kwargs
    ) -> np.ndarray:
        """
        %(base_model_ci.summary)s

        Parameters
        ----------
        %(base_model_ci.parameters)s
        level

        Returns
        -------
        %(base_model_ci.returns)s
        """  # noqa

        from rpy2 import robjects
        from rpy2.robjects import pandas2ri

        # TODO: DRY
        if self.model is None:
            raise RuntimeError(
                "Trying to call an uninitialized model. To initialize it, run `.fit()` first."
            )
        if self._lib is None:
            raise RuntimeError(
                f"Unable to run prediction, R package `{self._lib_name!r}` is not imported."
            )

        if not (0 <= level <= 1):
            raise ValueError(
                f"Expected level to be in interval `[0, 1]`, found `{level}`."
            )

        newdata = self._get_x_test(x_test)

        pandas2ri.activate()
        res = robjects.r.predict(
            self.model,
            newdata=pandas2ri.py2rpy(newdata),
            se=True,
        )
        pandas2ri.deactivate()

        self._y_test = np.array(res.rx2("fit")).squeeze().astype(self._dtype)
        se = np.array(res.rx2("se.fit")).squeeze().astype(self._dtype)

        level = norm.ppf(level + (1 - level) / 2)
        self._conf_int = np.c_[self.y_test - level * se, self.y_test + level * se]

        return self._conf_int

    @d.dedent
    def copy(self) -> "GAMR":
        """%(copy)s"""  # noqa
        res = GAMR(
            self.adata,
            self._n_splines,
            distribution=self._family,
            offset=self._original_offset,
            perform_import_check=False,
        )
        res._lib = self._lib
        res._lib_name = self._lib_name
        res._gam_kwargs = self._gam_kwargs

        return res

    def __getstate__(self) -> dict:
        return {k: v for k, v in self.__dict__.items() if k != "_lib"}

    def __setstate__(self, state: dict):
        self.__dict__ = state
        self._lib, self._lib_name = _maybe_import_r_lib(self._lib_name, raise_exc=True)


def _maybe_import_r_lib(
    name: str, raise_exc: bool = False
) -> Optional[Tuple[Any, str]]:
    global _r_lib, _r_lib_name

    if name == _r_lib_name and _r_lib is not None:
        return _r_lib, _r_lib_name

    try:
        from logging import ERROR

        from rpy2.robjects import r
        from rpy2.robjects.packages import PackageNotInstalledError, importr
        from rpy2.rinterface_lib.callbacks import logger

        logger.setLevel(ERROR)
        r["options"](warn=-1)

    except ImportError as e:
        raise ImportError(
            "Unable to import `rpy2`, install it first as `pip install rpy2` version `>=3.3.0`."
        ) from e

    try:
        _r_lib = importr(name)
        _r_lib_name = name
        return _r_lib, _r_lib_name
    except PackageNotInstalledError as e:
        if not raise_exc:
            return
        raise RuntimeError(
            f"Install R library `{name!r}` first as `install.packages({name!r}).`"
        ) from e
