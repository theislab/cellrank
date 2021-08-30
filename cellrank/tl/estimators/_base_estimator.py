from typing import Any, Dict, Tuple, Union, Mapping, Callable, Optional, Sequence
from typing_extensions import Literal

from abc import ABC, abstractmethod
from copy import copy
from copy import copy as copy_
from copy import deepcopy
from inspect import Parameter, signature, getmembers, currentframe
from contextlib import contextmanager

from anndata import AnnData
from cellrank._key import Key
from cellrank.ul._docs import d
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._lineage import Lineage
from cellrank.tl.estimators.mixins import IOMixin, KernelMixin, AnnDataMixin
from cellrank.tl.kernels._base_kernel import Kernel, KernelExpression

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix, csr_matrix


@d.get_sections(base="base_estimator", sections=["Parameters"])
class BaseEstimator(IOMixin, KernelMixin, AnnDataMixin, ABC):
    """
    Base class for all estimators.

    Parameters
    ----------
    obj
        Either a :class:`cellrank.tl.kernels.Kernel` object, an :class:`anndata.AnnData` object which
        stores the transition matrix in ``.obsp`` attribute or :mod:`numpy` or :mod:`scipy` array.
    obsp_key
        Key in :attr:`anndata.AnnData.obsp` where the transition matrix is stored.
        Only used when ``obj`` is an :class:`anndata.AnnData` object.
    kwargs
        Additional keyword arguments.
    """

    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        obsp_key: Optional[str] = None,
        **kwargs: Any,
    ):
        if isinstance(obj, Kernel):
            if obj._transition_matrix is None:
                raise ValueError(
                    "Compute transition matrix first as `.compute_transition_matrix()`."
                )
            kernel = obj
        elif isinstance(obj, KernelExpression):
            # this will fail if not all kernels have transition matrix computed
            kernel = obj.compute_transition_matrix()
        elif isinstance(obj, (np.ndarray, spmatrix)):
            kernel = PrecomputedKernel(obj)
        elif isinstance(obj, AnnData):
            if obsp_key is None:
                raise ValueError(
                    "Specify `obsp_key=...` when supplying an `AnnData` object."
                )
            elif obsp_key not in obj.obsp:
                raise KeyError(
                    f"Unable to find transition matrix in `adata.obsp[{obsp_key!r}]`."
                )
            kernel = PrecomputedKernel(obsp_key, adata=obj)
        else:
            raise TypeError(
                f"Expected an object of type `KernelExpression`, `numpy.ndarray`, `scipy.sparse.spmatrix` "
                f"or `anndata.AnnData`, got `{type(obj).__name__}`."
            )

        super().__init__(kernel=kernel, **kwargs)

        self._params: Dict[str, Any] = {}
        self._shadow_adata = AnnData(
            X=csr_matrix(self.adata.shape, dtype=self.adata.X.dtype),
            obs=self.adata.obs[[]].copy(),
            var=self.adata.var[[]].copy(),
            raw=None if self.adata.raw is None else self.adata.raw.to_adata(),
        )

    def __init_subclass__(cls, **kwargs: Any):
        super().__init_subclass__()

    def _set(
        self,
        attr: Optional[str] = None,
        obj: Optional[Union[pd.DataFrame, Mapping[str, Any]]] = None,
        key: Optional[str] = None,
        value: Optional[
            Union[np.ndarray, pd.Series, pd.DataFrame, Lineage, AnnData]
        ] = None,
        copy: bool = True,
        shadow_only: bool = False,
    ):
        """
        Set an attribute and optionally update ``obj[key]``.

        Parameters
        ----------
        attr
            Attribute to set. Only updated when we're not in the shadow. If `None`, don't update anything.
            See :attr:`_in_shadow` and ``obj`` for more information.
        obj
            Object which to update with ``value`` alongside the ``attr`.
            Usually, an attribute of :attr:`adata` is passed here.
        key
            Key in ``obj`` to update with ``value``. Only used when ``obj != None``.
        value
            Value to set. If `None` and ``key != None``, it removes the values under ``obj[key]``, if present.
        copy
            Whether to copy the ``value`` before setting it in ``obj``.
        shadow_only
            Whether or not to update the ``obj`` if we are not in the shadow.

        Returns
        -------
        Nothing, just optionally updates ``attr`` and/or ``obj[key]``.

        Raises
        ------
        AttributeError
            If ``attr`` doesn't exist.
        """
        if not self._in_shadow:
            if attr is not None:
                if not hasattr(self, attr):
                    raise AttributeError(attr)
                setattr(self, attr, value)
            if shadow_only:
                return

        if obj is None:
            return

        if key is not None:
            if value is None:
                try:
                    del obj[key]
                except KeyError:
                    pass
            else:
                obj[key] = copy_(value) if copy else value

    def _get(
        self,
        attr: str,
        obj: Union[pd.DataFrame, Mapping[str, Any]],
        key: str,
        where: Optional[Literal["obs", "obsm", "var", "varm", "uns"]] = None,
        dtype: Optional[Union[type, Tuple[type, ...]]] = None,
        copy: bool = True,
        allow_missing: bool = False,
    ) -> None:
        """
        Get data from an object and set an attribute.

        Parameters
        ----------
        attr
            Attribute to set.
        obj
            Object from which to extract the data.
        key
            Key in ``obj`` where the data is stored.
        where
            Attribute of :attr:`_shadow_adata` where to save the extracted data. If `None`, don't update it.
        dtype
            Valid type(s) of the extracted data.
        copy
            Copy the data before setting the ``attr``.
        allow_missing
            Whether or not to allow ``key`` to be missing in ``obj``.

        Returns
        -------
        Nothing, just updates ``attr`` with the extracted values and optionally :attr:`_shadow_adata`.

        Raises
        ------
        AttributeError
            If ``attr`` doesn't exist.
        TypeError
            If ``dtype != None`` and the extracted values are not instances of ``dtype``.
        KeyError
            If ``allow_missing = False`` and ``key`` was not found in ``obj``.
        """
        if not hasattr(self, attr):
            raise AttributeError(attr)

        try:
            data = obj[key]
            if dtype is not None and not isinstance(data, dtype):
                raise TypeError(
                    f"Expected `.{attr}` to be of type `{dtype}`, found `{type(data).__name__}`."
                )
            if copy:
                data = copy_(data)
            setattr(self, attr, data)
            if where is not None:
                getattr(self._shadow_adata, where)[key] = data
        except KeyError:
            if not allow_missing:
                raise

    @property
    @contextmanager
    def _shadow(self) -> None:
        """
        Temporarily set :attr:`adata` to :attr:`_shadow_adata`.

        Used to construct the serialization object in :meth:`to_adata`.
        """
        if self._in_shadow:
            yield
        else:
            adata = self.adata
            try:
                self.adata = self._shadow_adata
                yield
            finally:
                self.adata = adata

    @property
    def _in_shadow(self) -> bool:
        """Return `True` if :attr:`adata` is :attr:`_shadow_adata`."""
        return self.adata is self._shadow_adata

    def _create_params(
        self,
        locs: Optional[Mapping[str, Any]] = None,
        func: Optional[Callable] = None,
        remove: Sequence[str] = (),
    ) -> Dict[str, Any]:
        """
        Create parameters of interest from a function call.

        Parameters
        ----------
        locs
            Environment from which to get the parameters. If `None`, get the caller's environment.
        func
            Function of interest. If `None`, use the caller.
        remove
            Keys in ``locs`` which should not be included in the result.

        Returns
        -------
        The parameters as a :class:`dict`.

        Notes
        -----
        *args/**kwargs are always ignored and the values in ``locs`` are not copied.
        """
        frame = currentframe()
        try:
            if locs is None:
                locs = frame.f_back.f_locals
            if func is None or True:
                name = frame.f_back.f_code.co_name
                func = dict(getmembers(self)).get(name, None)
            if not callable(func):
                raise TypeError(
                    f"Expected `func` to be `callable`, found `{type(func).__name__}`."
                )

            params = {}
            for name, param in signature(func).parameters.items():
                if param.kind in (Parameter.VAR_POSITIONAL, Parameter.VAR_KEYWORD):
                    continue
                if name in remove:
                    continue
                if name in locs:
                    params[name] = locs[name]

            return params
        except AttributeError:
            # frame can be None
            return {}
        except TypeError:
            return {}
        finally:
            del frame

    def _read_params(self, key: str) -> Dict[str, Any]:
        """
        Read ``key`` from estimator params in :attr:`adata`.

        Usually called in :meth:`_read_adata` during :meth:`from_adata`.
        """
        ekey = Key.uns.estimator(self.backward) + "_params"
        return dict(self.adata.uns.get(ekey, {}).get(key, {}))

    # TODO: allow for keys to be copied from original?
    @d.dedent
    def to_adata(self) -> AnnData:
        """%(to_adata.full_desc)s"""  # noqa: D400
        adata = self._shadow_adata.copy()
        try:
            # kernel and estimator share the adata
            self.adata = adata
            self.kernel.write_to_adata()
        finally:
            self.adata = self.adata
        key = Key.uns.estimator(self.backward) + "_params"
        adata.uns[key] = deepcopy(self.params)
        return adata

    @classmethod
    @d.dedent
    def from_adata(cls, adata: AnnData, obsp_key: str) -> "BaseEstimator":
        """
        %(from_adata.full_desc)s

        Parameters
        ----------
        %(adata)s
        obsp_key
            Key in :attr:`anndata.AnnData.obsp` where the transition matrix is stored.

        Returns
        -------
        %(from_adata.returns)s
        """  # noqa: D400
        return super().from_adata(adata, obsp_key=obsp_key)

    @d.dedent
    def copy(self, *, deep: bool = False) -> "BaseEstimator":
        """
        Return a copy of self.

        Parameters
        ----------
        deep
            Whether to return a deep copy or not. If `True`, this also copies the :attr:`adata`.

        Returns
        -------
        A copy of self.
        """
        k = deepcopy(self.kernel) if deep else copy(self.kernel)
        res = type(self)(k)
        for k, v in self.__dict__.items():
            if isinstance(v, Mapping):
                res.__dict__[k] = deepcopy(v)
            elif k != "_kernel":
                res.__dict__[k] = deepcopy(v) if deep else copy(v)

        return res

    def __copy__(self) -> "BaseEstimator":
        return self.copy(deep=False)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}[n={len(self)}, kernel={repr(self.kernel)}]"

    def __str__(self) -> str:
        return f"{self.__class__.__name__}[n={len(self)}, kernel={str(self.kernel)}]"

    @property
    def params(self) -> Dict[str, Any]:
        """Estimator parameters."""
        return self._params

    @abstractmethod
    def fit(self, *args: Any, **kwargs: Any) -> "BaseEstimator":
        """
        Fit an estimator.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        Self.
        """

    @abstractmethod
    def predict(self, *args: Any, **kwargs: Any) -> None:
        """
        Run a prediction.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        Nothing.
        """
