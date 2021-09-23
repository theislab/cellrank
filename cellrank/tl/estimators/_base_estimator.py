from typing import Any, Dict, List, Tuple, Union, Mapping, Callable, Optional, Sequence
from typing_extensions import Literal

from abc import ABC, abstractmethod
from copy import copy as copy_
from copy import deepcopy
from inspect import Parameter, signature, getmembers, currentframe
from contextlib import contextmanager

from anndata import AnnData
from cellrank._key import Key
from cellrank.ul._docs import d
from cellrank.tl._mixins import IOMixin, KernelMixin, AnnDataMixin
from cellrank.tl.kernels import PrecomputedKernel
from cellrank.tl._lineage import Lineage
from cellrank.tl.kernels._base_kernel import Kernel, KernelExpression

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix, csr_matrix

Attr_t = (
    Literal["X", "raw", "layers", "obs", "var", "obsm", "varm", "obsp", "varp", "uns"],
)


@d.get_sections(base="base_estimator", sections=["Parameters"])
class BaseEstimator(IOMixin, KernelMixin, AnnDataMixin, ABC):
    """
    Base class for all estimators.

    Parameters
    ----------
    obj
        Can be one of the following:

            - :class:`cellrank.tl.kernels.Kernel` - kernel object.
            - :class:`anndata.AnnData` - annotated data object containing transition matrix in
              :attr:`anndata.AnnData.obsp`.
            - :class:`numpy.ndarray` - row-normalized sparse transition matrix.
            - :class:`scipy.sparse.spmatrix` - row-normalized sparse transition matrix.
    obsp_key
        Key in :attr:`anndata.AnnData.obsp` where the transition matrix is stored.
        Only used when ``obj`` is an :class:`anndata.AnnData` object.
    """

    def __init__(
        self,
        obj: Union[AnnData, np.ndarray, spmatrix, KernelExpression],
        obsp_key: Optional[str] = None,
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

        super().__init__(kernel=kernel)

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
            Union[np.ndarray, pd.Series, pd.DataFrame, Lineage, AnnData, Dict[str, Any]]
        ] = None,
        copy: bool = True,
        shadow_only: bool = False,
    ) -> None:
        """
        Set an attribute and optionally update ``obj[{key}]``.

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
        Nothing, just optionally updates ``attr`` and/or ``obj[{key}]``.

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

    @d.dedent
    def to_adata(
        self,
        keep: Union[Literal["all"], Sequence[Attr_t]] = ("X", "raw"),
        *,
        copy: Union[bool, Sequence[Attr_t]] = True,
    ) -> AnnData:
        """
        %(to_adata.full_desc)s

        Parameters
        ----------
        keep
            Which attributes to keep from the underlying :attr:`adata`. Valid options are:

                - `'all'` - keep all attributes specified in the signature.
                - :class:`typing.Sequence` - keep only subset of these attributes.
                - :class:`dict` - the keys correspond the attribute names and values to a subset of keys
                  which to keep from this attribute. If the values are specified either as `True` or `'all'`,
                  everything from this attribute will be kept.
        copy
            Whether to copy the data. Can be specified on per-attribute basis. Useful for attributes that store arrays.
            Attributes not specified here will not be copied.

        Returns
        -------
        %(adata)s
        """  # noqa: D400

        def handle_attribute(attr: Attr_t, keys: List[str], *, copy: bool) -> None:
            try:
                if attr == "X":
                    adata.X = deepcopy(self.adata.X) if copy else self.adata.X
                    return
                if attr == "raw":
                    adata.raw = (
                        self.adata.raw.to_adata()
                        if self.adata.raw is not None
                        else None
                    )
                    return

                old = getattr(self.adata, attr)
                new = getattr(adata, attr)
                if keys == ["all"]:
                    keys = list(old.keys())

                # fmt: off
                if isinstance(new, pd.DataFrame):
                    old = old[keys]
                    # avoid duplicates
                    old = old[old.columns.difference(new.columns)]
                    setattr(adata, attr, pd.merge(new, old, how="inner", left_index=True, right_index=True, copy=copy))
                elif isinstance(new, Mapping):
                    old = {k: old[k] for k in keys}
                    # old has preference, since it's user supplied
                    setattr(adata, attr, {**new, **(deepcopy(old) if copy else old)})
                else:
                    raise TypeError(f"Expected `adata.{attr}` to be either `Mapping` or `pandas. DataFrame`, "
                                    f"found `{type(new).__name__}`.")
                # fmt: on
            except KeyError:
                missing = sorted(k for k in keys if k not in old)
                raise KeyError(
                    f"Unable to find key(s) `{missing}` in `adata.{attr}`."
                ) from None

        adata = self._shadow_adata.copy()
        _adata = self.adata
        try:
            # kernel and estimator share the adata
            self.adata = adata
            self.kernel.write_to_adata()
        finally:
            self.adata = _adata
        key = Key.uns.estimator(self.backward) + "_params"
        adata.uns[key] = deepcopy(self.params)

        # fmt: off
        if isinstance(keep, str):
            if keep == 'all':
                keep = ["X", "raw", "layers", "obs", "var", "obsm", "varm", "obsp", "varp", "uns"]
            else:
                keep = [keep]
        if not isinstance(keep, Mapping):
            keep = {attr: True for attr in keep}
        if isinstance(copy, bool):
            copy = {attr: copy for attr in keep}
        elif isinstance(copy, str):
            copy = [copy]
        if not isinstance(copy, Mapping):
            copy = {attr: True for attr in copy}
        # fmt: on

        for attr, keys in keep.items():
            if keys is True:
                keys = ["all"]
            elif isinstance(keys, str):
                keys = [keys]
            if keys is False or not len(keys):
                continue
            handle_attribute(attr, keys=list(keys), copy=copy.get(attr, False))

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
        k = deepcopy(self.kernel) if deep else copy_(self.kernel)
        res = type(self)(k)
        for k, v in self.__dict__.items():
            if isinstance(v, Mapping):
                res.__dict__[k] = deepcopy(v)
            elif k != "_kernel":
                res.__dict__[k] = deepcopy(v) if deep else copy_(v)

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
