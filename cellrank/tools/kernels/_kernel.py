# -*- coding: utf-8 -*-
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
from functools import wraps, reduce

import numpy as np
from scipy.sparse import spdiags, issparse, spmatrix, csr_matrix

from scvelo.preprocessing.moments import get_moments

import jax.numpy as jnp
from jax import jacfwd, jacrev
from cellrank import logging as logg
from cellrank.tools._utils import (
    bias_knn,
    _normalize,
    _get_neighs,
    _has_neighs,
    _predict_fwd,
    _vals_to_csr,
    is_connected,
    is_symmetric,
)
from cellrank.utils._utils import _write_graph_data
from cellrank.tools._constants import Direction, _transition

_ERROR_DIRECTION_MSG = "Can only combine kernels that have the same direction."
_ERROR_EMPTY_CACHE_MSG = (
    "Fatal error: tried to used cached values, but the cache was empty."
)
_ERROR_CONF_ADAPT = (
    "Confidence adaptive operator is only supported for kernels, found type `{}`."
)
_ERROR_VAR_NOT_FOUND = "Variances not found in kernel `{}`."

_LOG_USING_CACHE = "Using cached transition matrix"

_n_dec = 2
_dtype = np.float64
_cond_num_tolerance = 1e-15
AnnData = TypeVar("AnnData")


class KernelExpression(ABC):
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

        If not present, compute it, if all the underlying kernels have been initialized.
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
    def adata(self) -> AnnData:
        """The annotated data object."""  # noqa
        pass

    @property
    def params(self) -> Dict[str, Any]:
        """Parameters which are used to compute the transition matrix."""
        return self._params

    def _format_params(self):
        return ", ".join(
            f"{k}={round(v, _n_dec) if isinstance(v, float) else v}"
            for k, v in self.params.items()
        )

    @transition_matrix.setter
    def transition_matrix(self, value: Union[np.ndarray, spmatrix]) -> None:
        """
        Set a new value of the transition matrix.

        Params
        ------
        value
            The new transition matrix. If the expression has no parent, the matrix is normalized.

        Returns
        -------
        None
        """

        if self._parent is None:
            self._transition_matrix = _normalize(value)
        else:
            self._transition_matrix = _normalize(value) if self._normalize else value

    @abstractmethod
    def compute_transition_matrix(self, *args, **kwargs) -> "KernelExpression":
        """
        Compute a transition matrix.

        Params
        ------
        args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        :class:`cellrank.tl.kernels.KernelExpression`
            Self.
        """
        pass

    def write_to_adata(self, key_added: Optional[str] = None) -> None:
        """
        Write the parameters and transition matrix to the underlying :paramref:`adata` object.

        Params
        ------
        key_added
            Postfix to be added to :paramref:`adata` `.uns` when writing the parameters.

        Returns
        -------
        None
            Updates the underlying :paramref:`.adata` object with the following fields:

                - .obsp[`'T_{fwd,bwd}` _ :paramref:`key_added` '] - transition matrix.
                - .uns[`'T_{fwd,bwd}` _  :paramref:`key_added` `_params'`] - parameters used.
        """

        if self.transition_matrix is None:
            raise ValueError(
                "Compute transition matrix first as `.compute_transition_matrix()`.`"
            )

        key = _transition(self._direction)
        if key_added is not None:
            key += f"_{key_added}"

        self.adata.uns[f"{key}_params"] = str(self)
        _write_graph_data(self.adata, self.transition_matrix, key)

    @abstractmethod
    def copy(self) -> "KernelExpression":
        """Return a copy of itself. Note that the underlying :paramref:`adata` object is not copied."""
        pass

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
        return list(self._get_kernels())

    def __xor__(self, other: "KernelExpression") -> "KernelExpression":
        return self.__rxor__(other)

    def __rxor__(self, other: "KernelExpression") -> "KernelExpression":
        def convert(obj):
            if isinstance(obj, _adaptive_add_type):
                if obj._variances is None:
                    raise ValueError(_ERROR_VAR_NOT_FOUND.format(obj))
                return KernelMul(
                    [
                        ConstantMatrix(
                            obj.adata, 1, obj._variances, backward=obj.backward
                        ),
                        obj,
                    ]
                )
            if isinstance(_is_bin_mult(obj, return_constant=False), _adaptive_add_type):
                e, c = _get_expr_and_constant(obj)
                if e._variances is None:
                    raise ValueError(_ERROR_VAR_NOT_FOUND.format(e))
                return KernelMul(
                    [ConstantMatrix(e.adata, c, e._variances, backward=e.backward), e]
                )

            return obj

        if (
            not isinstance(self, _adaptive_add_type)
            and not isinstance(self, KernelAdaptiveAdd)
            and not isinstance(
                _is_bin_mult(self, return_constant=False), _adaptive_add_type
            )
        ):
            raise TypeError(_ERROR_CONF_ADAPT.format(self.__class__.__name__))
        if (
            not isinstance(other, _adaptive_add_type)
            and not isinstance(other, KernelAdaptiveAdd)
            and not isinstance(
                _is_bin_mult(other, return_constant=False), _adaptive_add_type
            )
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


class UnaryKernelExpression(KernelExpression, ABC):
    """Base class for unary kernel expressions, such as kernels or constants."""

    def __init__(
        self,
        adata,
        backward: bool = False,
        op_name: Optional[str] = None,
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
        **kwargs,
    ):
        super().__init__(op_name, backward=backward, compute_cond_num=compute_cond_num)
        assert (
            op_name is None
        ), "Unary kernel does not support any kind operation associated with it."
        self._adata = adata
        self._read_from_adata(check_connectivity=check_connectivity, **kwargs)

    @abstractmethod
    def _read_from_adata(self, **kwargs):
        """Import the base-KNN graph and optionally check for symmetry and connectivity."""

        if not _has_neighs(self.adata):
            raise KeyError("Compute KNN graph first as `scanpy.pp.neighbors()`.")

        self._conn = _get_neighs(self.adata, "connectivities").astype(_dtype)

        check_connectivity = kwargs.pop("check_connectivity", False)
        if check_connectivity:
            start = logg.debug("Checking the KNN graph for connectedness")
            if not is_connected(self._conn):
                logg.warning("KNN graph is not connected", time=start)
            else:
                logg.debug("Knn graph is connected", time=start)

        start = logg.debug("Checking the KNN graph for symmetry")
        if not is_symmetric(self._conn):
            logg.warning("KNN graph is not symmetric", time=start)
        else:
            logg.debug("KNN graph is symmetric", time=start)

    def density_normalize(
        self, other: Union[np.ndarray, spmatrix]
    ) -> Union[np.ndarray, spmatrix]:
        """
        Density normalization by the underlying KNN graph.

        Params
        ------
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

    @property
    def adata(self) -> AnnData:
        """The annotated data object."""  # noqa
        return self._adata

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

    def _maybe_recalculate_constants(self, typp: Type):
        if typp == Constant:
            accessor = "transition_matrix"
        elif typp == ConstantMatrix:
            accessor = "_value"
        else:
            raise RuntimeError(
                f"Unable to determine accessor for type `{type.__name__}`."
            )

        constants = [_is_bin_mult(k, typp) for k in self]
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
    def adata(self) -> AnnData:
        """The annotated data object."""  # noqa
        # we can do this because Constant requires adata as well
        return self._kexprs[0].adata

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


class Kernel(UnaryKernelExpression, ABC):
    """
    A base class from which all kernels are derived.

    These kernels read from a given AnnData object, usually the KNN graph and additional variables, to compute a
    weighted, directed graph. Every kernel object has a direction. The kernels defined in the derived classes are not
    strictly kernels in the mathematical sense because they often only take one input argument - however, they build
    on other functions which have computed a similarity based on two input arguments. The role of the kernels defined
    here is to add directionality to these symmetric similarity relations or to transform them.

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    backward
        Direction of the process.
    kwargs
        Keyword arguments which can specify key to be read from :paramref:`adata` object.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        compute_cond_num: bool = False,
        **kwargs,
    ):
        super().__init__(
            adata, backward, op_name=None, compute_cond_num=compute_cond_num, **kwargs
        )

    def _get_kernels(self) -> Iterable["Kernel"]:
        yield self

    def _compute_transition_matrix(
        self, matrix: spmatrix, density_normalize: bool = True,
    ):

        # density correction based on node degrees in the KNN grpah
        if density_normalize:
            matrix = self.density_normalize(matrix)

        # check for zero-rows
        problematic_indices = np.where(np.array(matrix.sum(1)).flatten() == 0)[0]
        if len(problematic_indices) != 0:
            logg.warning(
                f"Detected {len(problematic_indices)} absorbing states in the transition matrix. "
                f"This matrix won't be reducible, consider setting `use_negative_cosines=True`"
            )
            for ix in problematic_indices:
                matrix[ix, ix] = 1.0

        # setting this property automatically row-normalizes
        self.transition_matrix = csr_matrix(matrix)
        self._maybe_compute_cond_num()


class Constant(Kernel):
    """Kernel representing a multiplication by a constant number."""

    def __init__(
        self, adata: AnnData, value: Union[int, float], backward: bool = False
    ):
        super().__init__(adata, backward=backward)
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"Value must be on `float` or `int`, found `{type(value).__name__}`."
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

    def copy(self) -> "Constant":
        """Return a copy of self."""
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
                f"Value must be on `float` or `int`, found `{type(value).__name__}`."
            )
        if value <= 0:
            raise ValueError(f"Expected the constant to be positive, found `{value}`.")

        self._value = value
        self._variances = (
            csr_matrix(variances) if not issparse(variances) else variances
        )
        self._recalculate(value)

    def copy(self) -> "ConstantMatrix":
        """Return a copy of self."""
        return ConstantMatrix(
            self.adata, self._value, copy(self._variances), self.backward
        )

    def _recalculate(self, value) -> None:
        self._value = value
        self._params = {"value": value}
        self._transition_matrix = value * self._variances

    def _read_from_adata(self, **kwargs) -> None:
        super()._read_from_adata(**kwargs)  # we need the shape info from neighbors

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


class PrecomputedKernel(Kernel):
    """Kernel which contains precomputed transition matrix."""

    def __init__(
        self,
        transition_matrix: Union[np.ndarray, spmatrix],
        adata: Optional[AnnData] = None,
        backward: bool = False,
        compute_cond_num: bool = False,
    ):
        from anndata import AnnData as _AnnData

        if not isinstance(transition_matrix, (np.ndarray, spmatrix)):
            raise TypeError(
                f"Expected transition matrix to be of type `numpy.ndarray` or `scipy.sparse.spmatrix`, "
                f"found `{type(transition_matrix).__name__}`."
            )

        if transition_matrix.shape[0] != transition_matrix.shape[1]:
            raise ValueError(
                f"Expected transition matrix to be square, found `{transition_matrix.shape}`."
            )

        if not np.allclose(np.sum(transition_matrix, axis=1), 1.0):
            raise ValueError("Not a valid transition matrix: not all rows sum to 1.")

        if adata is None:
            logg.debug("Creating empty dummy AnnData object")
            adata = _AnnData(
                csr_matrix((transition_matrix.shape[0], 1), dtype=np.float32)
            )

        super().__init__(adata, backward=backward, compute_cond_num=compute_cond_num)
        self._transition_matrix = csr_matrix(transition_matrix)
        self._maybe_compute_cond_num()

    def _read_from_adata(self, **kwargs):
        pass

    def copy(self) -> "PrecomputedKernel":
        """Return a copy of self."""
        pk = PrecomputedKernel(
            copy(self.transition_matrix), adata=self.adata, backward=self.backward
        )
        pk._cond_num = self.condition_number
        return pk

    def compute_transition_matrix(self, *args, **kwargs) -> "PrecomputedKernel":
        """Return self."""
        return self

    def __invert__(self) -> "PrecomputedKernel":
        # do not call parent's invert, since it removes the transition matrix
        self._direction = Direction.FORWARD if self.backward else Direction.BACKWARD
        return self

    def __repr__(self):
        return f"{'~' if self.backward and self._parent is None else ''}<Precomputed>"

    def __str__(self):
        return repr(self)


class VelocityKernel(Kernel):
    """
    Kernel which computes a transition matrix based on velocity correlations.

    This borrows ideas from both [Manno18]_ and [Bergen19]_. In short, for each cell *i*, we compute transition
    probabilities :math:`p_{i, j}` to each cell *j* in the neighborhood of *i*. The transition probabilities are
    computed as a multinominal logistic regression where the weights :math:`w_j` (for all *j*) are given by the vector
    that connects cell *i* with cell *j* in gene expression space, and the features :math:`x_i` are given
    by the velocity vector :math:`v_i` of cell *i*.

    Optionally, we apply a density correction as described in [Coifman05]_, where we use the implementation of
    [Haghverdi16]_.

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    backward
        Direction of the process.
    vkey
        Key in :paramref:`adata` `.uns` where the velocities are stored.
    var_key
        Key in :paramref:`adata` `.obsp` where the velocity variances are stored.
    use_negative_cosines
        Whether to use correlations with cells that have an angle > 90 degree with :math:`v_i`.
    compute_cond_num
        Whether to compute condition number of the transition matrix. Note that this might be costly,
        since it does not use sparse implementation.
    check_connectivity
        Check whether the underlying KNN graph is connected.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        vkey: str = "velocity",
        xkey: str = "Ms",
        gene_subset: Optional[Iterable] = None,
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            vkey=vkey,
            xkey=xkey,
            gene_subset=gene_subset,
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
        )
        self._vkey = vkey  # for copy
        self._xkey = xkey
        self._gene_subset = gene_subset

    def _read_from_adata(self, **kwargs):
        super()._read_from_adata(**kwargs)

        # check whether velocities have been computed
        vkey = kwargs.pop("vkey", "velocity")
        if vkey not in self.adata.layers.keys():
            raise KeyError("Compute RNA velocity first as `scv.tl.velocity()`.")

        # restrict genes to a subset, i.e. velocity genes or user provided list
        gene_subset = kwargs.pop("gene_subset", None)
        subset = np.ones(self.adata.n_vars, bool)
        if gene_subset is not None:
            var_names_subset = self.adata.var_names.isin(gene_subset)
            subset &= var_names_subset if len(var_names_subset) > 0 else gene_subset
        elif f"{vkey}_genes" in self.adata.var.keys():
            subset &= np.array(self.adata.var[f"{vkey}_genes"].values, dtype=bool)

        # chose data representation to use for transcriptomic displacements
        xkey = kwargs.pop("xkey", "Ms")
        xkey = xkey if xkey in self.adata.layers.keys() else "spliced"

        # filter both the velocities and the gene expression profiles to the gene subset. Densify the matrices.
        X = np.array(
            self.adata.layers[xkey].A[:, subset]
            if issparse(self.adata.layers[xkey])
            else self.adata.layers[xkey][:, subset]
        )
        V = np.array(
            self.adata.layers[vkey].A[:, subset]
            if issparse(self.adata.layers[vkey])
            else self.adata.layers[vkey][:, subset]
        )

        # remove genes that have any Nan values (in both X and V)
        nans = np.isnan(np.sum(V, axis=0))
        if np.any(nans):
            X = X[:, ~nans]
            V = V[:, ~nans]

        # check the velocity parameters
        if vkey + "_params" in self.adata.uns.keys():
            velocity_params = self.adata.uns[vkey + "_params"]
        else:
            velocity_params = None
            logg.debug("Unable to load velocity parameters. ")

        # add to self
        self._velocity = V
        self._gene_expression = X
        self._velocity_params = velocity_params

        logg.debug(
            "Adding `._velocity`, the gene-wise velocities. "
            "Adding `._gene_expression`, the gene expression values. "
            "Adding `._velocity_params`, the velocity parameters. "
        )

    def compute_transition_matrix(
        self,
        density_normalize: bool = False,
        backward_mode: str = "transpose",
        sigma_corr: float = 1.0,
        mode: str = "stochastic",
        use_negative_correlations: bool = True,
        **kwargs,
    ) -> "VelocityKernel":
        """
        Compute transition matrix based on velocity directions on the local manifold.

        For each cell, infer transition probabilities based on the correlation of the cell's
        velocity-extrapolated cell state with cell states of its *K* nearest neighbors.

        Params
        ------
        density_normalize
            Whether or not to use the underlying KNN graph for density normalization.
        backward_mode
            Options are `['transpose', 'negate']`. Only matters if initialized as :paramref:`backward` =`True`.
        sigma_corr
            Kernel width for exp kernel to be used to compute transition probabilities.
            from the velocity graph. If `None`, the median cosine correlation in absolute value is used.
        mode
            How to compute transition probabilities. Options are "stochastic" (propagate uncertainty analytically),
            "deterministic" (don't propagate uncertainty) and "sampling" (sample from velocity distribution)

        Returns
        -------
        None
            Nothing, but makes :paramref:`transition_matrix` available.
        """

        if mode not in ["stochastic", "deterministic", "sampling"]:
            raise NotImplementedError(f"Mode `{mode}` is not implemented. ")

        start = logg.info("Computing transition matrix based on velocity correlations")

        params = dict(  # noqa
            dnorm=density_normalize,
            bwd_mode=backward_mode if self._direction == Direction.BACKWARD else None,
            sigma_corr=sigma_corr,
            use_negative_correlations=use_negative_correlations,
            mode=mode,
        )

        # check whether we already computed such a transition matrix. If yes, load from cache
        if params == self._params:
            assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE, time=start)
            logg.info("    Finish", time=start)
            return self

        self._params = params

        # compute moments and define a function which returns hessian matrices
        if mode in ["stochastic", "sampling"]:
            velocity_expectation = get_moments(
                self.adata, self._velocity, second_order=False
            )
            velocity_variance = get_moments(
                self.adata, self._velocity, second_order=True
            )

        # define a function to compute hessian matrices
        get_hessian_fwd = jacfwd(jacrev(_predict_fwd))

        # loop over all cells
        vals, rows, cols, n_obs = [], [], [], self._gene_expression.shape[0]
        for i in range(10):
            print(f"i = {i}/{n_obs}")

            # get the neighbors
            nbhs_ixs = self._conn[i, :].indices
            W = self._gene_expression[nbhs_ixs, :] - self._gene_expression[nbhs_ixs, :]

            if mode == "deterministic":

                # evaluate the prediction at the actual velocity vector
                v_i = self._velocity[i, :]
                p = _predict_fwd(v_i, W, sigma=sigma_corr)

            elif mode == "stochastic":

                # get the expectation and variance
                v_i = velocity_expectation[i, :]
                variances = velocity_variance[i, :]

                # compute the Hessian tensor, and turn it into a matrix that has the diagonal elements in its rows
                H = get_hessian_fwd(v_i, W, sigma_corr)
                H_diag = jnp.concatenate([jnp.diag(h)[None, :] for h in H])

                # compute zero order term
                p_0 = _predict_fwd(v_i, W, sigma=sigma_corr)

                # compute second order term (note that the first order term cancels)
                p_2 = 0.5 * jnp.dot(H_diag, variances)

                # combine both to give the second order Taylor approximation
                p = p_0 + p_2

            elif mode == "sampling":
                pass

            # add the computed transition matrices to a list
            vals.extend(p)
            rows.extend(np.ones(len(nbhs_ixs)) * i)
            cols.extend(nbhs_ixs)

        # collect everything in a single transition matrix
        vals = np.hstack(vals)
        vals[np.isnan(vals)] = 0

        matrix = _vals_to_csr(vals, rows, cols, shape=(n_obs, n_obs))

        self._compute_transition_matrix(matrix, density_normalize=density_normalize)

        logg.info("    Finish", time=start)

        return self

    def copy(self) -> "VelocityKernel":
        """Return a copy of self."""
        vk = VelocityKernel(
            self.adata,
            backward=self.backward,
            vkey=self._vkey,
            xkey=self._xkey,
            gene_subset=self._gene_subset,
        )
        vk._params = copy(self.params)
        vk._cond_num = self.condition_number
        vk._transition_matrix = copy(self._transition_matrix)

        return vk


class ConnectivityKernel(Kernel):
    """
    Kernel which computes transition probabilities based on transcriptomic similarities.

    As a measure for transcriptomic similarity, we use the weighted KNN graph computed using
    :func:`scanpy.pp.neighbors`,see [Wolf18]_.
    By definition, the resulting transition matrix is symmetric and cannot be used to learn about the direction of the
    developmental process under consideration.
    However, the velocity-derived transition matrix can be combined with the similarity-based transition matrix as
    a means of regularization.

    Optionally, we apply a density correction as described in [Coifman05]_, where we use the implementation of
    [Haghverdi16]_.

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    backward
        Direction of the process.
    var_key
        Key in :paramref:`adata` `.uns` where the velocity variances are stored.
    compute_cond_num
        Whether to compute condition number of the transition matrix. Note that this might be costly,
        since it does not use sparse implementation.
    check_connectivity
        Check whether the underlying KNN graph is connected
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        var_key: Optional[str] = None,
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            var_key=var_key,
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
        )
        self._var_key = var_key

    def _read_from_adata(self, var_key: Optional[str] = None, **kwargs):
        super()._read_from_adata(var_key=var_key, **kwargs)

    def compute_transition_matrix(
        self, density_normalize: bool = True, **kwargs
    ) -> "ConnectivityKernel":
        """
        Compute transition matrix based on transcriptomic similarity.

        Uses symmetric, weighted KNN graph to compute symmetric transition matrix. The connectivities are computed
        using :func:`scanpy.pp.neighbors`. Depending on the parameters used there, they can be UMAP connectivities or
        gaussian-kernel-based connectivities with adaptive kernel width.

        Params
        ------
        density_normalize
            Whether or not to use the underlying KNN graph for density normalization.

        Returns
        -------
        None
            Nothing, but makes :paramref:`transition_matrix` available.
        """

        start = logg.info("Computing transition matrix based on connectivities")

        params = {"dnorm": density_normalize}
        if params == self.params:
            assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE)
            logg.info("    Finish", time=start)
            return self

        self._params = params
        self._compute_transition_matrix(
            matrix=self._conn.copy(), density_normalize=density_normalize
        )

        logg.info("    Finish", time=start)

        return self

    def copy(self) -> "ConnectivityKernel":
        """Return a copy of self."""
        ck = ConnectivityKernel(
            self.adata, backward=self.backward, var_key=self._var_key
        )
        ck._params = copy(self.params)
        ck._cond_num = self.condition_number
        ck._transition_matrix = copy(self._transition_matrix)

        return ck


class PalantirKernel(Kernel):
    """
    Kernel which computes transition probabilities in a similar way to *Palantir*, see [Setty19]_.

    *Palantir* computes a KNN graph in gene expression space and a pseudotime, which it then uses to direct the edges of
    the KNN graph, such that they are more likely to point into the direction of increasing pseudotime. To avoid
    disconnecting the graph, it does not remove all edges that point into the direction of decreasing pseudotime
    but keeps the ones that point to nodes inside a close radius. This radius is chosen according to the local density.

    The implementation presented here won't exactly reproduce the original *Palantir* algorithm (see below)
    but the results are qualitatively very similar.

    Optionally, we apply density correction as described in [Coifman05]_, where we use the implementation of
    [Haghverdi16]_.

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    backward
        Direction of the process.
    time_key
        Key in :paramref:`adata` `.obs` where the pseudotime is stored.
    var_key
        Key in :paramref:`adata` `.uns` where the velocity variances are stored.
    compute_cond_num
        Whether to compute condition number of the transition matrix. Note that this might be costly,
        since it does not use sparse implementation.
    """

    def __init__(
        self,
        adata: AnnData,
        backward: bool = False,
        time_key: str = "dpt_pseudotime",
        var_key: Optional[str] = None,
        compute_cond_num: bool = False,
        check_connectivity: bool = False,
    ):
        super().__init__(
            adata,
            backward=backward,
            time_key=time_key,
            var_key=var_key,
            compute_cond_num=compute_cond_num,
            check_connectivity=check_connectivity,
        )
        self._time_key = time_key
        self._var_key = var_key

    def _read_from_adata(self, var_key: Optional[str] = None, **kwargs):
        super()._read_from_adata(var_key=var_key, **kwargs)

        time_key = kwargs.pop("time_key", "dpt_pseudotime")
        if time_key not in self.adata.obs.keys():
            raise KeyError(f"Could not find time key in `adata.obs[{time_key!r}]`.")
        logg.debug("Adding `.pseudotime`")

        self._pseudotime = np.array(self.adata.obs[time_key]).astype(_dtype)

        if np.nanmin(self.pseudotime) < 0:
            raise ValueError(
                f"Minimum pseudotime must be non-negative, found {np.nanmin(self.pseudotime)}."
            )

    def compute_transition_matrix(
        self, k: int = 3, density_normalize: bool = True, **kwargs
    ) -> "PalantirKernel":
        """
        Compute transition matrix based on KNN graph and pseudotemporal ordering.

        This is a re-implementation of the Palantir algorithm by [Setty19]_.
        Note that this won't exactly reproduce the original Palantir results, for three reasons:

            - 1. Palantir computes the KNN graph in a scaled space of diffusion components.
            - 2. Palantir uses its own pseudotime to bias the KNN graph which is not implemented here.
            - 3. Palantir uses a slightly different mechanism to ensure the graph remains connected when removing edges that point into the "pseudotime past".

        If you would like to reproduce the original results, please use the original Palantir algorithm.

        Params
        ------
        k
            :paramref:`k` is the number of neighbors to keep for each node, regardless of pseudotime.
            This is done to ensure that the graph remains connected.
        density_normalize
            Whether or not to use the underlying KNN graph for density normalization.

        Returns
        -------
        None
            Nothing, but makes :paramref:`transition_matrix` available.
        """  # noqa

        start = logg.info("Computing transition matrix based on Palantir-like kernel")

        # get the connectivities and number of neighbors
        if (
            "neighbors" in self.adata.uns
            and "params" in self.adata.uns["neighbors"]
            and "n_neighbors" in self.adata.uns["neighbors"]["params"]
        ):
            n_neighbors = self.adata.uns["neighbors"]["params"]["n_neighbors"]
        else:
            logg.warning(
                "Could not find 'n_neighbors' in `adata.uns['neighbors']['params']`. Using an estimate"
            )
            n_neighbors = np.min(self._conn.sum(1))

        params = dict(k=k, dnorm=density_normalize, n_neighs=n_neighbors)  # noqa
        if params == self._params:
            assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE)
            logg.info("    Finish", time=start)
            return self

        self._params = params

        # handle backward case and run biasing function
        pseudotime = (
            np.max(self.pseudotime) - self.pseudotime
            if self._direction == Direction.BACKWARD
            else self.pseudotime
        )
        biased_conn = bias_knn(
            conn=self._conn, pseudotime=pseudotime, n_neighbors=n_neighbors, k=k
        ).astype(_dtype)
        # make sure the biased graph is still connected
        if not is_connected(biased_conn):
            logg.warning("Biased KNN graph is disconnected")

        self._compute_transition_matrix(
            matrix=biased_conn, density_normalize=density_normalize
        )

        logg.info("    Finish", time=start)

        return self

    @property
    def pseudotime(self) -> np.array:
        """Pseudotemporal ordering of cells."""
        return self._pseudotime

    def copy(self) -> "PalantirKernel":
        """Return a copy of self."""
        pk = PalantirKernel(
            self.adata,
            backward=self.backward,
            time_key=self._time_key,
            var_key=self._var_key,
        )
        pk._pseudotime = copy(self.pseudotime)
        pk._params = copy(self._params)
        pk._cond_num = self.condition_number
        pk._transition_matrix = copy(self._transition_matrix)

        return pk


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
            if kexpr.transition_matrix is None:
                if isinstance(kexpr, Kernel):
                    raise RuntimeError(
                        f"Kernel `{kexpr}` is uninitialized. "
                        f"Compute its transition matrix as `.compute_transition_matrix()`."
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

    def copy(self) -> "SimpleNaryExpression":
        """Return a copy of self."""
        sne = SimpleNaryExpression(
            [k.copy() for k in self], op_name=self._op_name, fn=self._fn
        )
        sne._transition_matrix = copy(self._transition_matrix)

        return sne


class KernelAdd(SimpleNaryExpression):
    """Base class that represents the addition of :class:`KernelExpression`."""

    def __init__(self, kexprs: List[KernelExpression], op_name: str):
        super().__init__(kexprs, op_name=op_name, fn=_reduce(np.add, 0))


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
        super().__init__(kexprs, op_name="*", fn=_reduce(np.multiply, 1))


def _reduce(func: Callable, initial: Union[int, float]) -> Callable:
    """
    Wrap :func:`reduce` function for a given function and an initial state.

    Params
    ------
    func
        Function to be used in the reduction.
    initial
        Initial value for the reduction.

    Returns
    -------
    Callable
        The wrapped :func:`reduce` function.
    """

    @wraps(func)
    def wrapper(seq: Iterable):
        return reduce(func, seq, initial)

    return wrapper


def _get_expr_and_constant(k: KernelMul) -> Tuple[KernelExpression, Union[int, float]]:
    """
    Get the value of a constant in binary multiplication.

    Params
    ------
    k
        Binary multiplication involving a constant and a kernel.

    Returns
    -------
    :class:`KernelExpression`, Union[int, float]
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

    Params
    ------
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
    :class:`cellrank.tl..kernels.KernelExpression`
        Depending on :paramref:`return_constant`, it either returns the constant multiplier
        or the expression being multiplied.
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


_adaptive_add_type = (
    ConnectivityKernel,
    VelocityKernel,
    PalantirKernel,
    PrecomputedKernel,
)
