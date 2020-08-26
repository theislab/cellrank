# -*- coding: utf-8 -*-
"""Module containing all models somehow interfacing R."""
from typing import Any, Tuple, Optional

import numpy as np
import pandas as pd

from cellrank.ul._docs import d
from cellrank.ul.models import BaseModel
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
    smoothing_param
        Smoothing parameter. Increasing this increases the smootheness of splines.
    distribution
        Family in `rpy2.robjects.r`, such as `'gaussian'` or `'poisson'`.
    backend
        R library used to fit GAMs. Valid options are `'mgcv'` or `'gam'`.
        Option `'gam'` ignores the number of splines, as well as ``distribution``and ``smoothing_param``.
    """  # noqa

    _fallback_backends = {
        "gam": "mgcv",
        "mgcv": "gam",
    }

    def __init__(
        self,
        adata: AnnData,
        n_splines: int = 5,
        smoothing_param: float = 2,
        distribution: str = "gaussian",
        backend: str = "mgcv",
        perform_import_check: bool = True,
    ):
        super().__init__(adata, model=None)
        self._n_splines = n_splines
        self._sp = smoothing_param
        self._lib = None
        self._lib_name = None
        self._family = distribution

        if backend not in self._fallback_backends.keys():
            raise ValueError(
                f"Invalid backend library `{backend!r}`. Valid options are `{list(self._fallback_backends.keys())}`."
            )

        if (
            perform_import_check
        ):  # it's a bit costly to import, copying just passes the reference
            self._lib, self._lib_name = _maybe_import_r_lib(
                backend
            ) or _maybe_import_r_lib(self._fallback_backends[backend], raise_exc=True)

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

        from rpy2 import robjects
        from rpy2.robjects import Formula, pandas2ri

        super().fit(x, y, w, **kwargs)

        use_ixs = self.w > 0
        self._x = self.x[use_ixs]
        self._y = self.y[use_ixs]
        self._w = self.w[use_ixs]

        family = getattr(robjects.r, self._family, None)
        if family is None:
            family = robjects.r.gaussian

        pandas2ri.activate()
        df = pandas2ri.py2rpy(pd.DataFrame(np.c_[self.x, self.y], columns=["x", "y"]))

        if self._lib_name == "mgcv":
            self._model = self._lib.gam(
                Formula(f'y ~ s(x, k={self._n_splines}, bs="cs")'),
                data=df,
                sp=self._sp,
                family=family,
                weights=pd.Series(self.w),
            )
        elif self._lib_name == "gam":
            self._model = self._lib.gam(
                Formula("y ~ s(x)"),
                data=df,
                sp=self._sp,
                family=family,
                weights=pd.Series(self.w),
            )
        else:
            raise NotImplementedError(
                f"No fitting implemented for R library `{self._lib_name!r}`."
            )

        pandas2ri.deactivate()

        return self

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
                f"Unable to fit the model, R package `{self._lib_name!r}` is not imported."
            )

        x_test = self._check(key_added, x_test)

        pandas2ri.activate()
        self._y_test = (
            np.array(
                robjects.r.predict(
                    self.model,
                    newdata=pandas2ri.py2rpy(pd.DataFrame(x_test, columns=["x"])),
                )
            )
            .squeeze()
            .astype(self._dtype)
        )
        pandas2ri.deactivate()

        return self.y_test

    @d.dedent
    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs
    ) -> np.ndarray:
        """
        %(base_model_ci.summary)s

        This method uses the :meth:`default_confidence_interval`.

        Parameters
        ----------
        %(base_model_ci.parameters)s

        Returns
        -------
        %(base_model_ci.returns)s
        """  # noqa
        return self.default_confidence_interval(x_test=x_test, **kwargs)

    @d.dedent
    def copy(self) -> "GAMR":
        """%(copy)s"""  # noqa
        res = GAMR(
            self.adata,
            self._n_splines,
            self._sp,
            distribution=self._family,
            perform_import_check=False,
        )
        res._lib = self._lib
        res._lib_name = self._lib_name
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
