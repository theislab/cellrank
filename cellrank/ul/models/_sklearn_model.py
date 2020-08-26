# -*- coding: utf-8 -*-
"""Module containing model which wraps around :mod:`sklearn` estimators."""
from copy import deepcopy
from typing import Iterable, Optional
from inspect import signature

import numpy as np
from sklearn.base import BaseEstimator

from cellrank.ul._docs import d
from cellrank.ul.models import BaseModel
from cellrank.ul.models._base_model import AnnData


@d.dedent
class SKLearnModel(BaseModel):
    """
    Wrapper around almost any :mod:`sklearn` model.

    Parameters
    ----------
    %(adata)s
    model
        Instance of :mod:`sklearn` model.
    weight_name
        Name of the weight argument for :paramref:`model` ``.fit``.
    ignore_raise
        Do not raise an exception if weight argument is not found in the fitting function of :paramref:`model`.
        This is useful in case when weight is passed in ``**kwargs`` and cannot be determined from signature.
    """

    _fit_names = ("fit", "__init__")
    _predict_names = ("predict", "__call__")
    _weight_names = ("w", "weights", "sample_weight", "sample_weights")
    _conf_int_names = ("conf_int", "confidence_intervals")

    def __init__(
        self,
        adata: AnnData,
        model: BaseEstimator,
        weight_name: Optional[str] = None,
        ignore_raise: bool = False,
    ):
        if not isinstance(model, BaseEstimator):
            raise TypeError(
                f"Expected model to be of type `BaseEstimator`, found `{type(model).__name__!r}`."
            )

        super().__init__(adata, model)

        fit_name = self._find_func(self._fit_names)
        predict_name = self._find_func(self._predict_names)
        ci_name = self._find_func(self._conf_int_names, use_default=True, default=None)
        self._weight_name = None

        if weight_name is None:
            self._weight_name = self._find_arg_name(fit_name, self._weight_names)
        else:
            params = signature(getattr(self.model, fit_name)).parameters
            if not ignore_raise and weight_name not in params:
                raise ValueError(
                    f"Unable to detect `{weight_name!r}` in the signature of `{fit_name!r}`."
                    f"If it's in `kwargs`, set `ignore_raise=True.`"
                )
            self._weight_name = weight_name

        if self._weight_name is None:
            raise RuntimeError(
                f"Unable to determine weights for function `{fit_name!r}`, searched `{self._weight_names}`. "
                f"Consider specifying it manually as `weight_name=...`."
            )

        self._fit_fn = getattr(self.model, fit_name)
        self._pred_fn = getattr(self.model, predict_name)
        self._ci_fn = None if ci_name is None else getattr(self.model, ci_name, None)

    @d.dedent
    def fit(
        self,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        w: Optional[np.ndarray] = None,
        **kwargs,
    ) -> "SKLearnModel":
        """
        %(base_model_fit.full_desc)s

        Parameters
        ----------
        %(base_model_fit.parameters)s

        Returns
        -------
        :class:`cellrank.ul.models.SKLearnModel`
            Fits the model and returns self.
        """  # noqa

        super().fit(x, y, w, **kwargs)

        if self._weight_name is not None:
            kwargs[self._weight_name] = self._w

        self._model = self._fit_fn(self.x, self.y, **kwargs)

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

        x_test = self._check(key_added, x_test)

        self._y_test = self._pred_fn(x_test, **kwargs)
        self._y_test = np.squeeze(self._y_test)

        return self.y_test

    @d.dedent
    def confidence_interval(
        self, x_test: Optional[np.ndarray] = None, **kwargs
    ) -> np.ndarray:
        """
        %(base_model_ci.full_desc)s

        Parameters
        ----------
        %(base_model_ci.parameters)s

        Returns
        -------
        %(base_model_ci.returns)s
        """  # noqa

        if self._ci_fn is None:
            return self.default_confidence_interval(x_test=x_test, **kwargs)

        x_test = self._check("_x_test", x_test)
        self._conf_int = self._ci_fn(x_test, **kwargs)

        return self.conf_int

    def _find_func(
        self,
        func_names: Iterable[str],
        use_default: bool = False,
        default: Optional[str] = None,
    ) -> Optional[str]:
        """
        Find a function in :parmref:`model` from given names.

        If `None` is found, use :parmaref;`default` or raise a :class:`RuntimeError`.

        Parameters
        ----------
        func_names
            Function names to search. The first one found is returned.
        use_default
            Whether to return the default is it's `None` or raise :class:`RuntimeError`.
        default
            The default functional name to use if `None` was found.

        Returns
        -------
        str, None
            Name of the function or the default name.
        """

        for name in func_names:
            if hasattr(self.model, name) and callable(getattr(self.model, name)):
                return name
        if use_default:
            return default
        raise RuntimeError(
            f"Unable to find function and no default specified, searched for `{list(func_names)}`."
        )

    def _find_arg_name(
        self, func_name: Optional[str], param_names: Iterable[str]
    ) -> Optional[str]:
        """
        Find an argument in :paramref:`model`'s ``func_name``.

        Parameters
        ----------
        func_name
            Function name of :paramref:`model`.
        param_names
            Parameter names to search. The first one found is returned.

        Returns
        -------
        str or None
            The parameter name or `None`, if `None` was found or ``func_name`` was `None`.
        """

        if func_name is None:
            return None

        for param in signature(getattr(self.model, func_name)).parameters:
            if param in param_names:
                return param

        return None

    @property
    def model(self) -> BaseEstimator:
        """The underlying :mod:`sklearn` model."""  # noqa
        return self._model

    @d.dedent
    def copy(self) -> "SKLearnModel":
        """%(copy)s"""  # noqa
        return SKLearnModel(self.adata, deepcopy(self._model))
