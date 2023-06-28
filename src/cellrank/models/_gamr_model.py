import copy
import enum
from typing import Any, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scipy.stats as st

from anndata import AnnData

from cellrank._utils._docs import d, inject_docs
from cellrank._utils._enum import ModeEnum
from cellrank.models import BaseModel, FailedModel
from cellrank.models._utils import _OFFSET_KEY, _get_knotlocs, _get_offset

__all__ = ["GAMR"]

# cache so we don't have to re-load
_r_lib = None
_r_lib_name = None


class KnotLocs(ModeEnum):
    AUTO = enum.auto()
    DENSITY = enum.auto()


@inject_docs(key=_OFFSET_KEY, kloc=KnotLocs)
@d.dedent
class GAMR(BaseModel):
    """Wrapper around R's `mgcv <https://cran.r-project.org/web/packages/mgcv/>`__ package for fitting GAMs.

    Parameters
    ----------
    %(adata)s
    n_knots
        Number of knots.
    distribution
        Distribution family in `rpy2.robjects.r`, such as ``'gaussian'`` or ``'nb'`` for negative binomial.
        If ``'nb'``, raw count data in :attr:`~anndata.AnnData.raw` is always used.
    basis
        Basis for the smoothing term.
        See `here <https://www.rdocumentation.org/packages/mgcv/versions/1.8-33/topics/s>`__ for valid options.
    knotlocs
        Position of the knots. Can be one of the following:

        - ``{kloc.AUTO!r}`` - let `mgcv` handle the knot positions.
        - ``{kloc.DENSITY!r}`` - position the knots based on the density of the pseudotime.
    offset
        Offset term for the GAM. Only available when ``distribution='nb'``. If `'default'`, it is calculated
        according to :cite:`robinson:10`. The values are saved in :attr:`adata.obs['{key}'] <anndata.AnnData.obs>`.
        If :obj:`None`, no offset is used.
    smoothing_penalty
        Penalty for the smoothing term. The larger the value, the smoother the fitted curve.
    kwargs
        Keyword arguments for `gam control
        <https://www.rdocumentation.org/packages/mgcv/versions/1.8-33/topics/gam.control>`__.
    """

    def __init__(
        self,
        adata: AnnData,
        n_knots: int = 5,
        distribution: str = "gaussian",
        basis: str = "cr",
        knotlocs: Literal["auto", "density"] = KnotLocs.AUTO,
        offset: Optional[Union[np.ndarray, Literal["default"]]] = "default",
        smoothing_penalty: float = 1.0,
        **kwargs: Any,
    ):
        if n_knots <= 0:
            raise ValueError(f"Expected `n_splines` to be positive, found `{n_knots}`.")
        if smoothing_penalty < 0:
            raise ValueError(f"Expected `smoothing_penalty` to be non-negative, found `{smoothing_penalty}`.")

        super().__init__(adata, model=None)
        self._n_knots = n_knots
        self._family = distribution

        self._formula = f"y ~ s(x, bs={basis!r}, k={self._n_knots}, sp={smoothing_penalty})"
        self._design_mat = None
        self._offset = None
        self._knotslocs = KnotLocs(knotlocs)

        self._control_kwargs = copy.copy(kwargs)

        self._lib = None
        self._lib_name = None

        # it's a bit costly to import, copying just passes the reference
        if kwargs.pop("perform_import_check", True):
            self._lib, self._lib_name = _maybe_import_r_lib("mgcv")

        if distribution == "nb" and offset is not None:
            if not isinstance(offset, (np.ndarray, str)):
                raise TypeError(
                    f"Expected `offset` to be either `'default'` or `numpy.ndarray`," f"got `{type(offset).__name__}`."
                )

            if isinstance(offset, str):
                if offset != "default":
                    raise ValueError("Only value `'default'` is allowed when `offset` is a string.")
                offset = _get_offset(adata, use_raw=True, recompute=False)

            offset = np.asarray(offset, dtype=self._dtype)

            if offset.shape != (adata.n_obs,):
                raise ValueError(f"Expected offset to be of shape `{(adata.n_obs,)}`, found `{offset.shape}`.")
            adata.obs[_OFFSET_KEY] = offset

            self._formula += " + offset(offset)"
            self._offset = offset

    @d.dedent
    def prepare(
        self,
        *args: Any,
        **kwargs: Any,
    ) -> "GAMR":
        """%(base_model_prepare.full_desc)s

        This also removes the zero and negative weights and prepares the design matrix.

        Parameters
        ----------
        %(base_model_prepare.parameters)s

        Returns
        -------
        %(base_model_prepare.returns)s
        """  # noqa

        if self._family == "nb":
            kwargs["use_raw"] = True
        prepared = super().prepare(*args, **kwargs)
        if isinstance(prepared, FailedModel):
            return prepared

        use_ixs = self.w > 0
        self._x = self.x[use_ixs]
        self._y = self.y[use_ixs]
        self._w = self.w[use_ixs]

        self._design_mat = pd.DataFrame(
            np.c_[self.x, self.y],
            columns=["x", "y"],
        )

        if self._offset is not None:
            mask = np.isin(self.adata.obs_names, self._obs_names[use_ixs])
            self._design_mat["offset"] = self._offset[mask]

        return self

    @d.dedent
    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs: Any,
    ) -> "GAMR":
        """%(base_model_fit.full_desc)s

        Parameters
        ----------
        %(base_model_fit.parameters)s

        Returns
        -------
        Fits the model and returns self. Updates the following fields by filtering out :math:`0` weights :attr:`w`:

        - :attr:`x` - %(base_model_x.summary)s
        - :attr:`y` - %(base_model_y.summary)s
        - :attr:`w` - %(base_model_w.summary)s
        """  # noqa

        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter

        super().fit(x, y, w, **kwargs)

        # order seems to matter, pandas2ri + default_converter results in:
        # Conversion 'py2rpy' not defined for objects of type '<class 'rpy2.rlike.container.OrdDict'>'
        with localconverter(pandas2ri.converter + ro.default_converter):
            family = getattr(ro.r, self._family)
            kwargs = {}
            if self._knotslocs != KnotLocs.AUTO:
                kwargs["knots"] = pd.DataFrame(
                    _get_knotlocs(
                        self.x,
                        self._n_knots,
                        uniform=False,
                    ),
                    columns=["x"],
                )

            self._model = self._lib.gam(
                ro.Formula(self._formula),
                data=self._design_mat,
                family=family,
                weights=pd.Series(self.w),
                control=self._lib.gam_control(**self._control_kwargs),
                **kwargs,
            )

        return self

    def _get_x_test(self, x_test: Optional[np.ndarray] = None, key_added: str = "_x_test") -> pd.DataFrame:
        newdata = pd.DataFrame(self._check(key_added, x_test), columns=["x"])

        if "offset" in self._design_mat:
            newdata["offset"] = np.mean(self._design_mat["offset"])

        return newdata

    @d.dedent
    def predict(
        self,
        x_test: Optional[np.ndarray] = None,
        key_added: str = "_x_test",
        level: Optional[float] = None,
        **kwargs,
    ) -> np.ndarray:
        """%(base_model_predict.full_desc)s

        This method can also compute the confidence interval.

        Parameters
        ----------
        %(base_model_predict.parameters)s
        level
            Confidence level for confidence interval calculation.
            If :obj:`None`, don't compute the confidence interval. Must be in :math:`[0, 1]`.

        Returns
        -------
        %(base_model_predict.returns)s
        """  # noqa
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter

        if self.model is None:
            raise RuntimeError("Trying to call an uninitialized model. To initialize it, run `.fit()` first.")
        if self._lib is None:
            raise RuntimeError(f"Unable to run prediction, R package `{self._lib_name!r}` is not imported.")

        if level is not None and not (0 <= level <= 1):
            raise ValueError(f"Expected the confidence level to be in interval `[0, 1]`, found `{level}`.")

        newdata = self._get_x_test(x_test)

        with localconverter(pandas2ri.converter + ro.default_converter):
            res = ro.r.predict(
                self.model,
                newdata=ro.pandas2ri.py2rpy(newdata),
                type="response",
                se=level is not None,
            )

        if level is None:
            self._y_test = np.array(res).squeeze().astype(self._dtype)
        else:
            self._y_test = np.array(res.rx2("fit")).squeeze().astype(self._dtype)
            se = np.array(res.rx2("se.fit")).squeeze().astype(self._dtype)

            level = st.norm.ppf(level + (1 - level) / 2)
            self._conf_int = np.c_[self.y_test - level * se, self.y_test + level * se]

        return self.y_test

    @d.dedent
    def confidence_interval(self, x_test: Optional[np.ndarray] = None, level: float = 0.95, **kwargs) -> np.ndarray:
        """%(base_model_ci.summary)s

        Internally, this method calls :meth:`~cellrank.models.GAMR.predict` to extract the confidence interval,
        if needed.

        Parameters
        ----------
        %(base_model_ci.parameters)s
        level
            Confidence level.

        Returns
        -------
        %(base_model_ci.returns)s
        """  # noqa

        # this is 2x as fast as opposed to calling `robjects.r.predict` again
        # on my PC (Michal):
        #   -`.predict` take ~6.5ms (without se=True, it's ~5.5ms, 20% slowdown)
        #   -`.predict` withouty se=True + `.confidence_interval` take ~12ms
        # this impl. achieves the 2x speedup withouty sacrificing the 20%

        if x_test is not None or level != 0.95 or self._conf_int is None:
            _ = self.predict(x_test=x_test, level=level)

        return self._conf_int

    def _deepcopy_attributes(self, dst: "GAMR") -> None:
        super()._deepcopy_attributes(dst)
        dst._offset = copy.deepcopy(self._offset)
        dst._design_mat = copy.deepcopy(self._design_mat)

    @d.dedent
    def copy(self) -> "GAMR":
        """%(copy)s"""  # noqa
        res = GAMR(
            self.adata,
            self._n_knots,
            distribution=self._family,
            offset=self._offset,
            knotlocs=self._knotslocs,
            perform_import_check=False,
        )
        self._shallowcopy_attributes(res)  # does not copy `prepared`, since we're not copying the arrays

        res._formula = self._formula  # we don't save the basis inside
        res._design_mat = self._design_mat
        res._offset = self._offset

        res._lib = self._lib  # just passing the reference
        res._lib_name = self._lib_name

        res._control_kwargs = self._control_kwargs

        return res

    def __getstate__(self) -> dict:
        return {k: v for k, v in self.__dict__.items() if k != "_lib"}

    def __setstate__(self, state: dict):
        self.__dict__ = state
        self._lib, self._lib_name = _maybe_import_r_lib(self._lib_name, raise_exc=True)


def _maybe_import_r_lib(name: str, raise_exc: bool = False) -> Tuple[Optional[Any], Optional[str]]:
    global _r_lib, _r_lib_name

    if name == _r_lib_name and _r_lib is not None:
        return _r_lib, _r_lib_name

    try:
        from logging import ERROR

        from rpy2.rinterface_lib.callbacks import logger
        from rpy2.robjects import r
        from rpy2.robjects.packages import PackageNotInstalledError, importr

        logger.setLevel(ERROR)
        r["options"](warn=-1)

    except ImportError as e:
        raise ImportError("Unable to import `rpy2`, install it first as `pip install rpy2` version `>=3.3.0`.") from e

    try:
        _r_lib = importr(name)
        _r_lib_name = name
        return _r_lib, _r_lib_name
    except PackageNotInstalledError as e:
        if not raise_exc:
            return None, None
        raise RuntimeError(f"Install R library `{name!r}` first as `install.packages({name!r}).`") from e
