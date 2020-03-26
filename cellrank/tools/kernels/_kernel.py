# -*- coding: utf-8 -*-
from abc import ABC, abstractmethod
from cellrank.tools._constants import Direction, _transition
from cellrank.tools._utils import (
    _normalize,
    is_connected,
    is_symmetric,
    bias_knn,
    has_neighs,
    get_neighs,
)

from typing import Optional, Union, Callable, List, Iterable, Tuple, Type, Any, Dict
from anndata import AnnData
from scanpy import logging as logg
from numpy import ndarray
from scipy.sparse import issparse, spdiags, csr_matrix, spmatrix
from functools import wraps, reduce
from copy import copy

import numpy as np

_ERROR_DIRECTION_MSG = "Can only combine kernels that have the same direction."
_ERROR_EMPTY_CACHE_MSG = (
    "Fatal error: tried to used cached values, but the cache was empty."
)
_ERROR_CONF_ADAPT = (
    "Confidence adaptive operator is only supported for kernels, found type `{}`."
)
_ERROR_VAR_NOT_FOUND = "Variances not found in kernel `{}`."

_LOG_USING_CACHE = "DEBUG: Using cached transition matrix"

_n_dec = 2
_dtype = np.float64


class KernelExpression(ABC):
    """
    Base class for all kernels and kernel expressions.
    """

    def __init__(self, op_name: Optional[str] = None, backward: bool = False):
        self._op_name = op_name
        self._transition_matrix = None
        self._direction = Direction.BACKWARD if backward else Direction.FORWARD
        self._params = dict()
        self._normalize = True
        self._parent = None

    @property
    def transition_matrix(self) -> Union[np.ndarray, spmatrix]:
        """
        Returns row-normalized transition matrix, if present, or tries computing it, if all underlying
        kernels have been initialized.
        """

        if self._parent is None and self._transition_matrix is None:
            self.compute_transition_matrix()

        return self._transition_matrix

    @property
    def backward(self) -> bool:
        """
        Return `True` if the direction of the process is backwards, otherwise `False`.
        """
        # to do proper error checking, we need to propagate this kind information
        return self._direction == Direction.BACKWARD

    @property
    @abstractmethod
    def adata(self) -> AnnData:
        """
        Get the annotated data object.

        Returns
        -------
        :class:`anndata.AnnData`
            The underlying :paramref:`.adata` object.
        """
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
            The new transition matrix.

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

    def write_to_adata(self, key_added: Optional[str] = None):
        """
        Write the parameters and transition matrix to the underlying adata object.

        Params
        ------
        key_added
            Postfix to be added to :paramref`.adata` `.uns.

        Returns
        -------
        None
            Updates the underlying :paramref:`.adata` object with the following:
                - `.uns[:paramref:`T_{fwd, bwd}` _`:paramref:`key_added`]['T']` - transition matrix
                - `.uns[:paramref:`T_{fwd, bwd}` _`:paramref:`key_added`]['params']` - parameters used for calculation
        """

        if self.transition_matrix is None:
            raise ValueError(
                "Compute transition matrix first as `.compute_transition_matrix()`.`"
            )

        key = _transition(self._direction)
        if key_added is not None:
            key += f"_{key_added}"

        if self.adata.uns.get(key, None) is not None:
            logg.debug(f"DEBUG: Overwriting key `{key!r}` in `adata.uns`")

        self.adata.uns[key] = dict()
        self.adata.uns[key]["params"] = str(self)
        self.adata.uns[key]["T"] = self.transition_matrix

        logg.debug(f"DEBUG: Added `{key!r}` to `adata.uns`")

    def __xor__(self, other) -> "KernelExpression":
        return self.__rxor__(other)

    def __rxor__(self, other) -> "KernelExpression":
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
            exprs = [k for k in s] + (
                [k for k in o] if isinstance(o, KernelAdaptiveAdd) else [o]
            )
            # (c1 * x + c2 * y + ...) + (c3 * z) => (c1 * x + c2 * y + c3 + ... + * z)
            if all((_is_bin_mult(k, ConstantMatrix) for k in exprs)):
                return KernelAdaptiveAdd(exprs)

        # same but reverse
        if isinstance(o, KernelAdaptiveAdd):
            exprs = ([k for k in s] if isinstance(s, KernelAdaptiveAdd) else [s]) + [
                k for k in o
            ]
            if all((_is_bin_mult(k, ConstantMatrix) for k in exprs)):
                return KernelAdaptiveAdd(exprs)

        return KernelAdaptiveAdd([s, o])

    def __add__(self, other) -> "KernelExpression":
        return self.__radd__(other)

    def __radd__(self, other) -> "KernelExpression":
        if not isinstance(other, KernelExpression):
            raise TypeError(
                f"Expected type `KernelExpression`, found `{other.__class__.__name__}`."
            )

        s = self * 1 if isinstance(self, Kernel) else self
        o = other * 1 if isinstance(other, Kernel) else other

        if isinstance(s, KernelSimpleAdd):
            exprs = [k for k in s] + [o]
            # (c1 * x + c2 * y + ...) + (c3 * z) => (c1 * x + c2 * y + c3 + ... + * z)
            if all((_is_bin_mult(k) for k in exprs)):
                return KernelSimpleAdd(exprs)

        # same but reverse
        if isinstance(o, KernelSimpleAdd):
            exprs = [s] + [k for k in o]
            if all((_is_bin_mult(k) for k in exprs)):
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

    def __rmul__(self, other) -> "KernelExpression":
        return self.__mul__(other)

    def __mul__(self, other) -> "KernelExpression":
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

    def __invert__(self) -> "KernelExpression":
        # mustn't return a copy because transition matrix
        self._transition_matrix = None
        self._params = dict()
        self._direction = Direction.FORWARD if self.backward else Direction.BACKWARD

        return self


class UnaryKernelExpression(KernelExpression, ABC):
    """
    Base class for unary kernel expressions, such as kernels or constants.
    """

    def __init__(
        self, adata, backward: bool = False, op_name: Optional[str] = None, **kwargs
    ):
        super().__init__(op_name, backward=backward)
        assert (
            op_name is None
        ), "Unary kernel does not support any kind operation associated with it."
        self._adata = adata
        self._read_from_adata(**kwargs)

    @abstractmethod
    def _read_from_adata(self, **kwargs):
        """
        Import the base-KNN graph and check for symmetry and connectivity.
        """

        if not has_neighs(self.adata):
            raise KeyError("Compute KNN graph first as `scanpy.pp.neighbors()`.")

        self._conn = get_neighs(self.adata, "connectivities").astype(_dtype)

        start = logg.debug("Checking the KNN graph for connectedness")
        if not is_connected(self._conn):
            logg.warning("KNN graph is not connected", time=start)

        start = logg.debug("Checking the KNN graph for symmetry")
        if not is_symmetric(self._conn):
            logg.warning("KNN graph is not symmetric", time=start)

        variance_key = kwargs.pop("variance_key", None)
        if variance_key is not None:
            logg.debug(f"DEBUG: Loading variances from `adata.uns[{variance_key!r}]`")
            variance_key = f"{variance_key}_variances"
            if variance_key in self.adata.uns.keys():
                # keep it sparse
                self._variances = csr_matrix(
                    self.adata.uns[variance_key].astype(_dtype)
                )
            else:
                self._variances = None
                logg.debug(
                    f"DEBUG: Unable to load variances`{variance_key}` from `adata.uns`"
                )
        else:
            logg.debug("DEBUG: No variance key specified")

    def density_normalize(
        self, other: Union[ndarray, spmatrix]
    ) -> Union[ndarray, spmatrix]:
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

        logg.debug("DEBUG: Density-normalizing the transition matrix")

        q = np.asarray(self._conn.sum(axis=0))

        if not issparse(other):
            Q = np.diag(1.0 / q)
        else:
            Q = spdiags(1.0 / q, 0, other.shape[0], other.shape[0])

        return Q @ other @ Q

    @property
    def adata(self):
        return self._adata

    def __repr__(self):
        return f"<{self.__class__.__name__[:4]}>"

    def __str__(self):
        return f"<{self.__class__.__name__[:4]}[{self._format_params()}]>"


class NaryKernelExpression(KernelExpression, ABC):
    """
    Base class for n-ary kernel expressions.
    """

    def __init__(self, kexprs: List[KernelExpression], op_name: Optional[str] = None):
        assert len(kexprs), "No kernel expressions specified."
        backward = kexprs[0].backward
        assert all((k.backward == backward for k in kexprs)), _ERROR_DIRECTION_MSG
        super().__init__(op_name, backward=backward)

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
        if all((c is not None for c in constants)):
            assert all(
                (isinstance(getattr(c, accessor), (int, float)) for c in constants)
            )
            total = sum((getattr(c, accessor) for c in constants)) + 0.0
            for c in constants:
                c._recalculate(getattr(c, accessor) / total)

            for kexpr in self._kexprs:  # don't normalize  (c * x)
                kexpr._normalize = False

    @property
    def adata(self):
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
            "("
            + f" {self._op_name} ".join((repr(kexpr) for kexpr in self._kexprs))
            + ")"
        )

    def __str__(self) -> str:
        return (
            "("
            + f" {self._op_name} ".join((str(kexpr) for kexpr in self._kexprs))
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

    def __init__(self, adata: AnnData, backward: bool = False, **kwargs):
        super().__init__(adata, backward, op_name=None, **kwargs)


class Constant(Kernel):
    """
    Class representing a multiplication by a constant number.
    """

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
        self._params = dict(value=value)

    def _read_from_adata(self, **kwargs):
        pass

    def compute_transition_matrix(self, *args, **kwargs) -> "Constant":
        return self

    def __invert__(self) -> "Constant":
        # do not call parent's invert, since it removes the transition matrix
        self._direction = Direction.FORWARD if self.backward else Direction.BACKWARD
        return self

    def __repr__(self) -> str:
        return repr(round(self.transition_matrix, _n_dec))

    def __str__(self) -> str:
        return str(round(self.transition_matrix, _n_dec))


class ConstantMatrix(Kernel):
    """
    Class representing multiplication by a constant matrix.
    """

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

    def _recalculate(self, value):
        self._value = value
        self._params = dict(value=value)
        self._transition_matrix = value * self._variances

    def _read_from_adata(self, **kwargs):
        super()._read_from_adata(**kwargs)  # we need the shape info from neighbors

    def compute_transition_matrix(self, *args, **kwargs) -> "KernelExpression":
        return self

    def __invert__(self) -> "ConstantMatrix":
        # do not call parent's invert, since it removes the transition matrix
        self._direction = Direction.FORWARD if self.backward else Direction.BACKWARD
        return self

    def __repr__(self) -> str:
        return repr(round(self._value, _n_dec))

    def __str__(self) -> str:
        return str(round(self._value, _n_dec))


class VelocityKernel(Kernel):
    """
    Implements a kernel class which computes a transition matrix based on velocity correlations.

    This borrows ideas from both [Manno18]_ and [Bergen19]_. In short, for each cell *i*, we compute transition
    probabilities :math:`p_{i, j}` to each cell *j* ( in the neighborhood of *i*. The transition probabilities are
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
    """

    def __init__(self, adata: AnnData, backward: bool = False, vkey: str = "velocity"):
        super().__init__(adata, backward=backward, vkey=vkey)

    def _read_from_adata(self, vkey: str, **kwargs):
        super()._read_from_adata(variance_key="velocity", **kwargs)
        if (vkey + "_graph" not in self.adata.uns.keys()) or (
            vkey + "_graph_neg" not in self.adata.uns.keys()
        ):
            raise ValueError(
                "Compute cosine correlations first as `scvelo.tl.velocity_graph()`."
            )

        velo_corr_pos, velo_corr_neg = (
            csr_matrix(self.adata.uns[vkey + "_graph"]).copy(),
            csr_matrix(self.adata.uns[vkey + "_graph_neg"]).copy(),
        )
        logg.debug("Adding `.velo_corr`, the velocity correlations")

        self.velo_corr = (velo_corr_pos + velo_corr_neg).astype(_dtype)

    def compute_transition_matrix(
        self,
        density_normalize: bool = True,
        backward_mode: str = "transpose",
        sigma_corr: Optional[float] = None,
        **kwargs,
    ) -> "VelocityKernel":
        """
        Compute transition matrix based on velocity correlations.

        For each cell, infer transition probabilities based on the correlation of the cell's
        velocity-extrapolated cell state with cell states of its K nearest neighbors.

        Params
        ------
        density_normalize
            Whether or not to use the underlying KNN graph for density normalization.
        backward_mode
            Options are `['transpose', 'negate']`. Only matters if initialised as :paramref:`backward` =`True`.
        sigma_corr
            Kernel width for exp kernel to be used to compute transition probabilities
            from the velocity graph. If `None`, the median cosine correlation in absolute value is used.

        Returns
        -------
        None
            Makes :paramref:`transition_matrix` available.
        """

        start = logg.info("Computing transition matrix based on velocity correlations")

        # get the correlations, handle backwards case
        if self._direction == Direction.BACKWARD:
            if backward_mode == "negate":
                correlations = self.velo_corr.multiply(-1)
            elif backward_mode == "transpose":
                correlations = self.velo_corr.T
            else:
                raise ValueError(f"Unknown backward mode `{backward_mode!r}`.")
        else:
            correlations = self.velo_corr

        # set the scaling parameter for the softmax
        med_corr = np.median(np.abs(correlations.data))
        if sigma_corr is None:
            sigma_corr = 1 / med_corr

        params = dict(
            dnorm=density_normalize,
            bwd_mode=backward_mode if self._direction == Direction.BACKWARD else None,
            sigma_corr=sigma_corr,
        )

        if params == self._params:
            assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE, time=start)
            logg.info("     Finish", time=start)
            return self

        self._params = params

        # compute directed graph --> multi class log reg
        velo_graph = correlations.copy()
        velo_graph.data = np.exp(velo_graph.data * sigma_corr)

        # normalize
        if density_normalize:
            velo_graph = self.density_normalize(velo_graph)
        logg.info("    Finish", time=start)

        self.transition_matrix = csr_matrix(velo_graph)

        return self


class ConnectivityKernel(Kernel):
    """
    Implements a kernel class which computes transition probabilities based on transcriptomic similarities.

    As a measure for transcriptomic similarity, we use the weighted KNN graph computed using :func:`scanpy.pp.neighbors`,
    see [Wolf18]_. By definition, the resulting transition matrix is symmetric and cannot be used to learn about
    the direction of the developmental process under consideration. However, the velocity-derived transition matrix can
    be combined with the similarity-based transition matrix as a means of regularization.

    Optionally, we apply a density correction as described in [Coifman05]_, where we use the implementation of
    [Haghverdi16]_.

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    backward
        Direction of the process.
    """

    def __init__(self, adata: AnnData, backward: bool = False):
        super().__init__(adata, backward=backward)

    def _read_from_adata(self, **kwargs):
        super()._read_from_adata(variance_key="connectivity", **kwargs)

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
            Makes :paramref:`transition_matrix` available.
        """

        start = logg.info("Computing transition matrix based on connectivities")

        params = dict(dnorm=density_normalize)
        if params == self._params:
            assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
            logg.debug(_LOG_USING_CACHE)
            logg.info("    Finish", time=start)
            return self

        self._params = params
        conn = self._conn.copy()

        if density_normalize:
            conn = self.density_normalize(conn)
        logg.info("    Finish", time=start)

        self.transition_matrix = csr_matrix(conn)

        return self


class PalantirKernel(Kernel):
    """
    Implements a kernel class which computes transition probabilities in a similar way to *Palantir*, see [Setty19]_

    *Palantir* computes a KNN graph in gene expression space and a pseudotime, which it then uses to direct the edges of
    the KNN graph, such that they are more likely to point into the direction of increasing pseudotime. To avoid
    disconnecting the graph, it does not remove all edges that point into the direction of decreasing pseudotime
    but keeps the ones that point to nodes inside a close radius. This radius is chosen according to the local density.

    The implementation presented here won't exactly reproduce the original *Palantir* algorithm (see below)
    but the results are qualitatively very similar.

    Optionally, we apply a density correction as described in [Coifman05]_,where we use the implementation of
    [Haghverdi16]_.

    Params
    ------
    adata : :class:`anndata.AnnData`
        Annotated data object.
    backward
        Direction of the process.
    time_key
        Key in :paramref:`adata` `.obs` where the pseudotime is stored.
    """

    def __init__(
        self, adata: AnnData, backward: bool = False, time_key: str = "dpt_pseudotime"
    ):
        super().__init__(adata, backward=backward, time_key=time_key)

    def _read_from_adata(self, time_key: str, **kwargs):
        super()._read_from_adata(variance_key="palantir", **kwargs)
        if time_key not in self.adata.obs.keys():
            raise KeyError(f"Could not find time key `{time_key!r}` in `adata.obs`.")
        logg.debug("Adding `.pseudotime`")

        self.pseudotime = np.array(self.adata.obs[time_key]).astype(_dtype)

        if np.min(self.pseudotime) < 0:
            raise ValueError(f"Pseudotime must be positive")

    def compute_transition_matrix(
        self, k: int = 3, density_normalize: bool = True, **kwargs
    ) -> "PalantirKernel":
        """
        Compute transition matrix based on KNN graph and pseudotemporal ordering.

        This is a re-implementation of the Palantir algorithm by [Setty19]_.
        Note that this won't exactly reproduce the original Palantir results, for three reasons:

        - 1. Palantir computes the KNN graph in a scaled space of diffusion components.
        - 2. Palantir uses its own pseudotime to bias the KNN graph which is not implemented here
        - 3. Palantir uses a slightly different mechanism to ensure the graph remains connected when removing edges that
             point into the "pseudotime past".
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
            Makes :paramref:`transition_matrix` available.
        """

        start = logg.info("Computing transition matrix based on Palantir-like kernel")

        # get the connectivities and number of neighbors
        if (
            "neighbors" in self.adata.uns.keys()
            and "params" in self.adata.uns["neighbors"]
            and "n_neighbors" in self.adata.uns["neighbors"]["params"].keys()
        ):
            n_neighbors = self.adata.uns["neighbors"]["params"]["n_neighbors"]
        else:
            logg.warning(
                f"Could not find 'n_neighbors' in `adata.uns['neighbors']['params']`. Using an estimate"
            )
            n_neighbors = np.min(self._conn.sum(1))

        params = dict(k=k, dnorm=density_normalize, n_neighs=n_neighbors)
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

        # normalize
        if density_normalize:
            biased_conn = self.density_normalize(biased_conn)
        logg.info("    Finish", time=start)

        self.transition_matrix = csr_matrix(biased_conn)

        return self


class SimpleNaryExpression(NaryKernelExpression):
    """
    Base class for n-ary operations.
    """

    def __init__(self, kexprs: List[KernelExpression], op_name: str, fn: Callable):
        super().__init__(kexprs, op_name=op_name)
        self._fn = fn

    def compute_transition_matrix(self, *args, **kwargs) -> "SimpleNaryExpression":
        # must be done before, because the underlying expression dont' have to be normed
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

        return self


class KernelAdd(SimpleNaryExpression):
    def __init__(self, kexprs: List[KernelExpression], op_name: str):
        super().__init__(kexprs, op_name=op_name, fn=_reduce(np.add, 0))


class KernelSimpleAdd(KernelAdd):
    """
    Addition between two kernel expressions.
    """

    def __init__(self, kexprs: List[KernelExpression]):
        super().__init__(kexprs, op_name="+")


class KernelAdaptiveAdd(KernelAdd):
    """
    Adaptive addition between two kernels.
    """

    def __init__(self, kexprs: List[KernelExpression]):
        super().__init__(kexprs, op_name="^")


class KernelMul(SimpleNaryExpression):
    """
    Multiplication between two kernel expressions.
    """

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
    :class:`KernelExpression` or Union[int, float]
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


_adaptive_add_type = (ConnectivityKernel, VelocityKernel, PalantirKernel)
