"""Kernel module."""

from abc import ABC, abstractmethod
from copy import copy
from typing import (
    Any,
    Dict,
    List,
    Type,
    Tuple,
    Union,
    TypeVar,
    Callable,
    Iterable,
    Optional,
)
from functools import reduce

import numpy as np
from scipy.sparse import spdiags, issparse, spmatrix, csr_matrix, isspmatrix_csr

from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._utils import (
    _connected,
    _normalize,
    _symmetric,
    _get_neighs,
    _has_neighs,
)
from cellrank.ul._utils import Pickleable, _write_graph_data
from cellrank.tl._constants import Direction, _transition
from cellrank.tl.kernels._utils import _filter_kwargs

_ERROR_DIRECTION_MSG = "Can only combine kernels that have the same direction."
_ERROR_EMPTY_CACHE_MSG = (
    "Fatal error: tried to used cached values, but the cache was empty."
)
_ERROR_CONF_ADAPT = (
    "Confidence adaptive operator is only supported for kernels, found type `{!r}`."
)
_ERROR_VAR_NOT_FOUND = "Variances not found in kernel `{!r}`."

_LOG_USING_CACHE = "Using cached transition matrix"

_RTOL = 1e-12
_n_dec = 2
_dtype = np.float64
_cond_num_tolerance = 1e-15
AnnData = TypeVar("AnnData")


class KernelExpression(Pickleable, ABC):
    """Base class for all kernels and kernel expressions."""

    def __init__(
        self,
        op_name: Optional[str] = None,
        backward: bool = False,
        compute_cond_num: bool = False,
    ):
        self._op_name = op_name
        self._transition_matrix = None
        self._direction = Direction.BACKWARD if backward else Direction.FORWARD
        self._compute_cond_num = compute_cond_num
        self._cond_num = None
        self._params = {}
        self._normalize = True
        self._parent = None

    @property
    def condition_number(self):
        """Condition number of the transition matrix."""
        return self._cond_num

    @property
    def transition_matrix(self) -> Union[np.ndarray, spmatrix]:
        """
        Return row-normalized transition matrix.

        If not present, it is computed, if all the underlying kernels have been initialized.
        """

        if self._parent is None and self._transition_matrix is None:
            self.compute_transition_matrix()

        return self._transition_matrix

    @property
    def backward(self) -> bool:
        """Direction of the process."""
        return self._direction == Direction.BACKWARD

    @property
    @abstractmethod
    @d.dedent
    def adata(self) -> AnnData:
        """
        Annotated data object.

        Returns
        -------
        %(adata_ret)s
        """

    @adata.setter
    @abstractmethod
    def adata(self, value):
        pass

    @property
    def params(self) -> Dict[str, Any]:
        """Parameters which are used to compute the transition matrix."""
        if len(self.kernels) == 1:
            return self._params
        # we need some identifier
        return {f"{repr(k)}:{i}": k.params for i, k in enumerate(self.kernels)}

    def _format_params(self):
        return ", ".join(
            f"{k}={round(v, _n_dec) if isinstance(v, float) else v}"
            for k, v in self.params.items()
        )

    @transition_matrix.setter
    def transition_matrix(self, value: Union[np.ndarray, spmatrix]) -> None:
        """
        Set a new value of the transition matrix.

        Parameters
        ----------
        value
            The new transition matrix. If the expression has no parent, the matrix is normalized, if needed.

        Returns
        -------
        None
            Nothing, just updates the :paramref:`transition_matrix` and optionally normalizes it.
        """
        should_norm = ~np.isclose(value.sum(1), 1.0, rtol=_RTOL).all()

        if self._parent is None:
            self._transition_matrix = _normalize(value) if should_norm else value
        else:
            # it's AND, not OR, because of combinations
            self._transition_matrix = (
                _normalize(value) if self._normalize and should_norm else value
            )

    @abstractmethod
    def compute_transition_matrix(self, *args, **kwargs) -> "KernelExpression":
        """
        Compute a transition matrix.

        Parameters
        ----------
        *args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        :class:`cellrank.tl.kernels.KernelExpression`
            Self.
        """

    @d.get_sections(base="write_to_adata", sections=["Parameters"])
    @inject_docs()  # get rid of {{}}
    @d.dedent
    def write_to_adata(self, key: Optional[str] = None) -> None:
        """
        Write the transition matrix and parameters used for computation to the underlying :paramref:`adata` object.

        Parameters
        ----------
        key
            Key used when writing transition matrix to :paramref:`adata`.
            If `None`, the ``key`` is set to `'T_bwd'` if :paramref:`backward` is `True`, else `'T_fwd'`.

        Returns
        -------
        None
            %(write_to_adata)s
        """

        if self._transition_matrix is None:
            raise ValueError(
                "Compute transition matrix first as `.compute_transition_matrix()`."
            )

        if key is None:
            key = _transition(self._direction)

        self.adata.uns[f"{key}_params"] = str(self)
        _write_graph_data(self.adata, self.transition_matrix, key)

    @abstractmethod
    def copy(self) -> "KernelExpression":
        """Return a copy of itself. Note that the underlying :paramref:`adata` object is not copied."""

    def _maybe_compute_cond_num(self):
        if self._compute_cond_num and self._cond_num is None:
            logg.debug(f"Computing condition number of `{repr(self)}`")
            self._cond_num = np.linalg.cond(
                self._transition_matrix.toarray()
                if issparse(self._transition_matrix)
                else self._transition_matrix
            )
            if self._cond_num > _cond_num_tolerance:
                logg.warning(
                    f"`{repr(self)}` may be ill-conditioned, its condition number is `{self._cond_num:.2e}`"
                )
            else:
                logg.info(f"Condition number is `{self._cond_num:.2e}`")

    @abstractmethod
    def _get_kernels(self) -> Iterable["Kernel"]:
        pass

    @property
    def kernels(self) -> List["Kernel"]:
        """Get the kernels of the kernel expression, except for constants."""
        return list(set(self._get_kernels()))

    def __xor__(self, other: "KernelExpression") -> "KernelExpression":
        return self.__rxor__(other)

    def __rxor__(self, other: "KernelExpression") -> "KernelExpression":
        def convert(obj):
            if _is_adaptive_type(obj):
                if obj._mat_scaler is None:
                    raise ValueError(_ERROR_VAR_NOT_FOUND.format(obj))
                return KernelMul(
                    [
                        ConstantMatrix(
                            obj.adata, 1, obj._mat_scaler, backward=obj.backward
                        ),
                        obj,
                    ]
                )
            if _is_adaptive_type(_is_bin_mult(obj, return_constant=False)):
                e, c = _get_expr_and_constant(obj)
                if e._mat_scaler is None:
                    raise ValueError(_ERROR_VAR_NOT_FOUND.format(e))
                return KernelMul(
                    [ConstantMatrix(e.adata, c, e._mat_scaler, backward=e.backward), e]
                )

            return obj

        if (
            not _is_adaptive_type(self)
            and not isinstance(self, KernelAdaptiveAdd)
            and not _is_adaptive_type(_is_bin_mult(self, return_constant=False))
        ):
            raise TypeError(_ERROR_CONF_ADAPT.format(self.__class__.__name__))
        if (
            not _is_adaptive_type(other)
            and not isinstance(other, KernelAdaptiveAdd)
            and not _is_adaptive_type(_is_bin_mult(other, return_constant=False))
        ):
            raise TypeError(_ERROR_CONF_ADAPT.format(other.__class__.__name__))

        s = convert(self)
        o = convert(other)

        if isinstance(s, KernelAdaptiveAdd):
            exprs = list(s) + (list(o) if isinstance(o, KernelAdaptiveAdd) else [o])
            # (c1 * x + c2 * y + ...) + (c3 * z) => (c1 * x + c2 * y + c3 + ... + * z)
            if all(_is_bin_mult(k, ConstantMatrix) for k in exprs):
                return KernelAdaptiveAdd(exprs)

        # same but reverse
        if isinstance(o, KernelAdaptiveAdd):
            exprs = (list(s) if isinstance(s, KernelAdaptiveAdd) else [s]) + list(o)
            if all(_is_bin_mult(k, ConstantMatrix) for k in exprs):
                return KernelAdaptiveAdd(exprs)

        return KernelAdaptiveAdd([s, o])

    def __add__(self, other: "KernelExpression") -> "KernelExpression":
        return self.__radd__(other)

    def __radd__(self, other: "KernelExpression") -> "KernelExpression":
        if not isinstance(other, KernelExpression):
            raise TypeError(
                f"Expected type `KernelExpression`, found `{other.__class__.__name__}`."
            )

        s = self * 1 if isinstance(self, Kernel) else self
        o = other * 1 if isinstance(other, Kernel) else other

        if isinstance(s, KernelSimpleAdd):
            exprs = list(s) + [o]
            # (c1 * x + c2 * y + ...) + (c3 * z) => (c1 * x + c2 * y + c3 + ... + * z)
            if all(_is_bin_mult(k) for k in exprs):
                return KernelSimpleAdd(exprs)

        # same but reverse
        if isinstance(o, KernelSimpleAdd):
            exprs = [s] + list(o)
            if all(_is_bin_mult(k) for k in exprs):
                return KernelSimpleAdd(exprs)

        # (c1 + c2) => c3
        if isinstance(s, Constant) and isinstance(o, Constant):
            assert s.backward == o.backward, _ERROR_DIRECTION_MSG
            return Constant(
                s.adata, s.transition_matrix + o.transition_matrix, backward=s.backward
            )

        ss = s * 1 if not isinstance(s, KernelMul) else s
        oo = o * 1 if not isinstance(o, KernelMul) else o

        return KernelSimpleAdd([ss, oo])

    def __rmul__(
        self, other: Union[int, float, "KernelExpression"]
    ) -> "KernelExpression":
        return self.__mul__(other)

    def __mul__(
        self, other: Union[float, int, "KernelExpression"]
    ) -> "KernelExpression":

        if isinstance(other, (int, float)):
            other = Constant(self.adata, other, backward=self.backward)

        if not isinstance(other, KernelExpression):
            raise TypeError(
                f"Expected type `KernelExpression`, found `{other.__class__.__name__}`."
            )

        # (c1 * c2) => c3
        if isinstance(self, Constant) and isinstance(
            other, Constant
        ):  # small optimization
            assert self.backward == other.backward, _ERROR_DIRECTION_MSG
            return Constant(
                self.adata,
                self.transition_matrix * other.transition_matrix,
                backward=self.backward,
            )

        s = (
            self
            if isinstance(self, (KernelMul, Constant))
            else KernelMul([Constant(self.adata, 1, backward=self.backward), self])
        )
        o = (
            other
            if isinstance(other, (KernelMul, Constant))
            else KernelMul([Constant(other.adata, 1, backward=other.backward), other])
        )

        cs, co = _is_bin_mult(s), _is_bin_mult(o)
        if cs and isinstance(o, Constant):
            cs._transition_matrix *= o.transition_matrix
            return s

        if co and isinstance(s, Constant):
            co._transition_matrix *= s.transition_matrix
            return o

        return KernelMul([s, o])

    def __invert__(self: "KernelExpression") -> "KernelExpression":
        # mustn't return a copy because transition matrix
        self._transition_matrix = None
        self._params = {}
        self._direction = Direction.FORWARD if self.backward else Direction.BACKWARD

        return self

    def __copy__(self) -> "KernelExpression":
        return self.copy()

    def __deepcopy__(self, memodict={}) -> "KernelExpression":  # noqa
        res = self.copy()
        if self._parent is None:
            res.adata = self.adata.copy()
        memodict[id(self)] = res

        return res


class UnaryKernelExpression(KernelExpression, ABC):
    """Base class for unary kernel expressions, such as kernels or constants."""

    def __init__(
        self,
        adata,
        backward: bool = False,
        op_name: Optional[str] = None,
        compute_cond_num: bool = False,
        **kwargs,
    ):
        super().__init__(op_name, backward=backward, compute_cond_num=compute_cond_num)
        assert (
            op_name is None
        ), "Unary kernel does not support any kind operation associated with it."
        self._adata = adata

    @property
    @d.dedent
    def adata(self) -> AnnData:
        """
        Annotated data object.

        Returns
        -------
        %(adata_ret)s
        """
        return self._adata

    @adata.setter
    def adata(self, _adata: AnnData):
        from anndata import AnnData

        if not isinstance(_adata, AnnData):
            raise TypeError(
                f"Expected argument of type `anndata.AnnData`, found `{type(_adata).__name__!r}`."
            )
        # otherwise, we'd have to reread bunch of attributes - it's better to initialize new object
        if _adata.shape != self.adata.shape:
            raise ValueError(
                f"Expected the new object to have same shape as previous object `{self.adata.shape}`, "
                f"found `{_adata.shape}`."
            )
        self._adata = _adata

    def __repr__(self):
        return f"{'~' if self.backward and self._parent is None else ''}<{self.__class__.__name__[:4]}>"

    def __str__(self):
        params_fmt = self._format_params()
        if params_fmt:
            return (
                f"{'~' if self.backward and self._parent is None else ''}"
                f"<{self.__class__.__name__[:4]}[{params_fmt}]>"
            )
        return repr(self)


class NaryKernelExpression(KernelExpression, ABC):
    """Base class for n-ary kernel expressions."""

    def __init__(self, kexprs: List[KernelExpression], op_name: Optional[str] = None):
        assert len(kexprs), "No kernel expressions specified."

        backward = kexprs[0].backward
        assert all(k.backward == backward for k in kexprs), _ERROR_DIRECTION_MSG

        # use OR instead of AND
        super().__init__(
            op_name,
            backward=backward,
            compute_cond_num=any(k._compute_cond_num for k in kexprs),
        )

        # copies of constants are necessary because of the recalculation
        self._kexprs = [copy(k) if isinstance(k, Constant) else k for k in kexprs]

        for kexprs in self._kexprs:
            kexprs._parent = self

    def _maybe_recalculate_constants(self, type_: Type):
        if type_ == Constant:
            accessor = "transition_matrix"
        elif type_ == ConstantMatrix:
            accessor = "_value"
        else:
            raise RuntimeError(
                f"Unable to determine accessor for type `{type.__name__!r}`."
            )

        constants = [_is_bin_mult(k, type_) for k in self]
        if all(c is not None for c in constants):
            assert all(
                isinstance(getattr(c, accessor), (int, float)) for c in constants
            )
            total = sum(getattr(c, accessor) for c in constants) + 0.0
            for c in constants:
                c._recalculate(getattr(c, accessor) / total)

            for kexpr in self._kexprs:  # don't normalize  (c * x)
                kexpr._normalize = False

    def _get_kernels(self) -> Iterable["Kernel"]:
        for k in self:
            if isinstance(k, Kernel) and not isinstance(k, (Constant, ConstantMatrix)):
                yield k
            elif isinstance(k, NaryKernelExpression):
                yield from k._get_kernels()

    @property
    @d.dedent
    def adata(self) -> AnnData:
        """
        Annotated data object.

        Returns
        -------
        %(adata_ret)s
        """
        # we can do this because Constant requires adata as well
        return self._kexprs[0].adata

    @adata.setter
    def adata(self, _adata: AnnData):
        self._kexprs[0].adata = _adata

    def __invert__(self) -> "NaryKernelExpression":
        super().__invert__()
        self._kexprs = [~kexpr for kexpr in self]
        return self

    def __len__(self) -> int:
        return len(self._kexprs)

    def __getitem__(self, item) -> "KernelExpression":
        return self._kexprs[item]

    def __repr__(self) -> str:
        return (
            f"{'~' if self.backward and self._parent is None else ''}("
            + f" {self._op_name} ".join(repr(kexpr) for kexpr in self._kexprs)
            + ")"
        )

    def __str__(self) -> str:
        return (
            f"{'~' if self.backward and self._parent is None else ''}("
            + f" {self._op_name} ".join(str(kexpr) for kexpr in self._kexprs)
            + ")"
        )


@d.dedent
class Kernel(UnaryKernelExpression, ABC):
    """
    A base class from which all kernels are derived.

    These kernels read from a given AnnData object, usually the KNN graph and additional variables, to compute a
    weighted, directed graph. Every kernel object has a direction. The kernels defined in the derived classes are not
    strictly kernels in the mathematical sense because they often only take one input argument - however, they build
    on other functions which have computed a similarity based on two input arguments. The role of the kernels defined
    here is to add directionality to these symmetric similarity relations or to transform them.

    Parameters
    ----------
    %(adata)s
    %(backward)s
    compute_cond_num
        Whether to compute the condition number of the transition matrix. For large matrices, this can be very slow.
    check_connectivity
        Check whether the underlying KNN graph is connected.
    kwargs
        Keyword arguments which can specify key to be read from :paramref:`adata` object.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
        **kwargs,
    ):
        super().__init__(
            adata, backward, op_name=None, compute_cond_num=compute_cond_num, **kwargs
        )
        self._read_from_adata(check_connectivity=check_connectivity, **kwargs)

    def _read_from_adata(self, **kwargs):
        """Import the base-KNN graph and optionally check for symmetry and connectivity."""

        if not _has_neighs(self.adata):
            raise KeyError("Compute KNN graph first as `scanpy.pp.neighbors()`.")

        self._conn = _get_neighs(self.adata, "connectivities").astype(_dtype)

        check_connectivity = kwargs.pop("check_connectivity", False)
        if check_connectivity:
            start = logg.debug("Checking the KNN graph for connectedness")
            if not _connected(self._conn):
                logg.warning("KNN graph is not connected", time=start)
            else:
                logg.debug("KNN graph is connected", time=start)

        start = logg.debug("Checking the KNN graph for symmetry")
        if not _symmetric(self._conn):
            logg.warning("KNN graph is not symmetric", time=start)
        else:
            logg.debug("KNN graph is symmetric", time=start)

    def _density_normalize(
        self, other: Union[np.ndarray, spmatrix]
    ) -> Union[np.ndarray, spmatrix]:
        """
        Density normalization by the underlying KNN graph.

        Parameters
        ----------
        other:
            Matrix to normalize.

        Returns
        -------
        :class:`np.ndarray` or :class:`scipy.sparse.spmatrix`
            Density normalized transition matrix.
        """

        logg.debug("Density-normalizing the transition matrix")

        q = np.asarray(self._conn.sum(axis=0))

        if not issparse(other):
            Q = np.diag(1.0 / q)
        else:
            Q = spdiags(1.0 / q, 0, other.shape[0], other.shape[0])

        return Q @ other @ Q

    def _get_kernels(self) -> Iterable["Kernel"]:
        yield self

    def _compute_transition_matrix(
        self,
        matrix: spmatrix,
        density_normalize: bool = True,
    ):
        # density correction based on node degrees in the KNN graph
        matrix = csr_matrix(matrix) if not isspmatrix_csr(matrix) else matrix

        if density_normalize:
            matrix = self._density_normalize(matrix)

        # check for zero-rows
        problematic_indices = np.where(np.array(matrix.sum(1)).flatten() == 0)[0]
        if len(problematic_indices):
            logg.warning(
                f"Detected `{len(problematic_indices)}` absorbing states in the transition matrix. "
                f"This matrix won't be reducible"
            )
            matrix[problematic_indices, problematic_indices] = 1.0

        # setting this property automatically row-normalizes
        self.transition_matrix = matrix
        self._maybe_compute_cond_num()


@d.dedent
class Constant(Kernel):
    """
    Kernel representing a multiplication by a constant number.

    Parameters
    ----------
    %(adata)s
    value
        Constant value by which to multiply Must be a positive number.
    %(backward)s
    """

    def __init__(
        self, adata: AnnData, value: Union[int, float], backward: bool = False
    ):
        super().__init__(adata, backward=backward)
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"Value must be a `float` or `int`, found `{type(value).__name__}`."
            )
        if value <= 0:
            raise ValueError(f"Expected the constant to be positive, found `{value}`.")
        self._recalculate(value)

    def _recalculate(self, value):
        self._transition_matrix = value
        self._params = {"value": value}

    def _read_from_adata(self, **kwargs):
        pass

    def compute_transition_matrix(self, *args, **kwargs) -> "Constant":
        """Return self."""
        return self

    @d.dedent
    def copy(self) -> "Constant":
        """%(copy)s"""  # noqa: D400, D401
        return Constant(self.adata, self.transition_matrix, self.backward)

    def __invert__(self) -> "Constant":
        # do not call parent's invert, since it removes the transition matrix
        self._direction = Direction.FORWARD if self.backward else Direction.BACKWARD
        return self

    def __repr__(self) -> str:
        return repr(round(self.transition_matrix, _n_dec))

    def __str__(self) -> str:
        return str(round(self.transition_matrix, _n_dec))


class ConstantMatrix(Kernel):
    """Kernel representing multiplication by a constant matrix."""

    def __init__(
        self,
        adata: AnnData,
        value: Union[int, float],
        variances: Union[np.ndarray, spmatrix],
        backward: bool = False,
    ):
        super().__init__(adata, backward=backward)
        conn_shape = self._conn.shape
        if variances.shape != conn_shape:
            raise ValueError(
                f"Expected variances of shape `{conn_shape}`, found `{variances.shape}`."
            )
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"Value must be on `float` or `int`, found `{type(value).__name__!r}`."
            )
        if value <= 0:
            raise ValueError(f"Expected the constant to be positive, found `{value}`.")

        self._value = value
        self._mat_scaler = (
            csr_matrix(variances) if not issparse(variances) else variances
        )
        self._recalculate(value)

    @d.dedent
    def copy(self) -> "ConstantMatrix":
        """%(copy)s"""  # noqa: D400, D401
        return ConstantMatrix(
            self.adata, self._value, copy(self._mat_scaler), self.backward
        )

    def _recalculate(self, value) -> None:
        self._value = value
        self._params = {"value": value}
        self._transition_matrix = value * self._mat_scaler

    def compute_transition_matrix(self, *args, **kwargs) -> "ConstantMatrix":
        """Return self."""
        return self

    def __invert__(self) -> "ConstantMatrix":
        # do not call parent's invert, since it removes the transition matrix
        self._direction = Direction.FORWARD if self.backward else Direction.BACKWARD
        return self

    def __repr__(self) -> str:
        return repr(round(self._value, _n_dec))

    def __str__(self) -> str:
        return str(round(self._value, _n_dec))


class SimpleNaryExpression(NaryKernelExpression):
    """Base class for n-ary operations."""

    def __init__(self, kexprs: List[KernelExpression], op_name: str, fn: Callable):
        super().__init__(kexprs, op_name=op_name)
        self._fn = fn

    def compute_transition_matrix(self, *args, **kwargs) -> "SimpleNaryExpression":
        """Compute and combine the transition matrices."""
        # must be done before, because the underlying expression don't have to be normed
        if isinstance(self, KernelSimpleAdd):
            self._maybe_recalculate_constants(Constant)
        elif isinstance(self, KernelAdaptiveAdd):
            self._maybe_recalculate_constants(ConstantMatrix)

        for kexpr in self:
            if kexpr._transition_matrix is None:
                if isinstance(kexpr, Kernel):
                    raise RuntimeError(
                        f"Kernel `{kexpr}` is uninitialized. "
                        f"Compute its transition matrix first as `.compute_transition_matrix()`."
                    )
                kexpr.compute_transition_matrix()
            elif isinstance(kexpr, Kernel):
                logg.debug(_LOG_USING_CACHE)

        self.transition_matrix = csr_matrix(
            self._fn([kexpr.transition_matrix for kexpr in self])
        )

        # only the top level expression and kernels will have condition number computed
        if self._parent is None:
            self._maybe_compute_cond_num()

        return self

    @d.dedent
    def copy(self) -> "SimpleNaryExpression":
        """%(copy)s"""  # noqa: D400, D401
        constructor = type(self)
        kwargs = {"op_name": self._op_name, "fn": self._fn}

        # preserve the type information so that combination can properly work
        # we test for this
        sne = constructor(
            [copy(k) for k in self], **_filter_kwargs(constructor, **kwargs)
        )

        # TODO: copying's a bit buggy - the parent stays the same
        # we could disallow it for inner expressions
        sne._transition_matrix = copy(self._transition_matrix)
        sne._condition_number = self._cond_num
        sne._normalize = self._normalize

        return sne


class KernelAdd(SimpleNaryExpression):
    """Base class that represents the addition of :class:`KernelExpression`."""

    def __init__(self, kexprs: List[KernelExpression], op_name: str):
        super().__init__(kexprs, op_name=op_name, fn=Reductor(np.add, 0))


class KernelSimpleAdd(KernelAdd):
    """Addition of :class:`KernelExpression`."""

    def __init__(self, kexprs: List[KernelExpression]):
        super().__init__(kexprs, op_name="+")


class KernelAdaptiveAdd(KernelAdd):
    """Adaptive addition of :class:`KernelExpression`."""

    def __init__(self, kexprs: List[KernelExpression]):
        super().__init__(kexprs, op_name="^")


class KernelMul(SimpleNaryExpression):
    """Multiplication of :class:`KernelExpression`."""

    def __init__(self, kexprs: List[KernelExpression]):
        super().__init__(kexprs, op_name="*", fn=Reductor(np.multiply, 1))


class Reductor:
    """
    Class that reduces an iterable.

    We don't define a decorator because of pickling.

    Parameters
    ----------
    func
        Reduction function.
    initial
        Initial value for :func:`reduce`.
    """

    def __init__(self, func: Callable, initial: Union[int, float]):
        self._func = func
        self._initial = initial

    def __call__(self, seq: Iterable):  # noqa
        return reduce(self._func, seq, self._initial)


def _get_expr_and_constant(k: KernelMul) -> Tuple[KernelExpression, Union[int, float]]:
    """
    Get the value of a constant in binary multiplication.

    Parameters
    ----------
    k
        Binary multiplication involving a constant and a kernel.

    Returns
    -------
    :class:`KernelExpression`, int or float
        The expression which is being multiplied and the value of the constant.
    """

    if not isinstance(k, KernelMul):
        raise TypeError(
            f"Expected expression to be of type `KernelMul`, found `{type(k).__name__}`."
        )
    if len(k) != 2:
        raise ValueError(
            f"Expected expression to be binary, found, `{len(k)}` subexpressions."
        )
    e1, e2 = k[0], k[1]

    if isinstance(e1, Constant):
        return e2, e1.transition_matrix
    elif isinstance(e2, Constant):
        return e1, e2.transition_matrix
    else:
        raise ValueError(
            "Expected one of the subexpressions to be `Constant`, found "
            f"`{type(e1).__name__}` and `{type(e2).__name__}`."
        )


def _is_bin_mult(
    k: KernelExpression,
    const_type: Union[Type, Tuple[Type]] = Constant,
    return_constant: bool = True,
) -> Optional[KernelExpression]:
    """
    Check if an expression is a binary multiplication.

    Parameters
    ----------
    k
        Kernel expression to check.
    const_type
        Type of the constant.
    return_constant
        Whether to return the constant in the expression or the multiplied expression.

    Returns
    -------
    None
        If the expression is not a binary multiplication.
    :class:`cellrank.tl.kernels.KernelExpression`
        Depending on ``return_constant``, it either returns the constant multiplier or the expression being multiplied.
    """

    if not isinstance(k, KernelMul):
        return None

    if len(k) != 2:
        return None

    lhs, rhs = k[0], k[1]

    if isinstance(lhs, const_type):
        return lhs if return_constant else rhs
    if isinstance(rhs, const_type):
        return rhs if return_constant else lhs

    return None


def _is_adaptive_type(k: KernelExpression) -> bool:
    return isinstance(k, Kernel) and not isinstance(k, (Constant, ConstantMatrix))
