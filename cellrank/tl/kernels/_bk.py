from typing import Any, Dict, List, Tuple, Union, Optional, Sequence

from abc import ABC, abstractmethod

from anndata import AnnData
from cellrank import logging as logg
from cellrank.ul._docs import d, inject_docs
from cellrank.tl._utils import _normalize
from cellrank.tl.kernels._mixins import BidirectionalMixin, UnidirectionalMixin

import numpy as np
from scipy.sparse import issparse, spmatrix, csr_matrix, isspmatrix_csr

Tmat_t = Union[np.ndarray, spmatrix]


class KernelExpression(ABC):
    def __init__(
        self,
        parent: Optional["KernelExpression"] = None,
        **_: Any,
    ):
        self._parent = parent
        self._normalize = parent is None
        self._transition_matrix = None
        self._params: Dict[str, Any] = {}

    @abstractmethod
    def compute_transition_matrix(
        self, *args: Any, **kwargs: Any
    ) -> "KernelExpression":
        """
        Compute a transition matrix.

        Parameters
        ----------
        args
            Positional arguments.
        kwargs
            Keyword arguments.

        Returns
        -------
        Modifies :attr:`transition_matrix` and returns self.
        """

    @abstractmethod
    def copy(self, *, deep: bool = False) -> "KernelExpression":
        """Return a copy of itself. The underlying :attr:`adata` object is not copied."""

    @property
    @abstractmethod
    def adata(self) -> AnnData:
        """Annotated data object."""

    @adata.setter
    @abstractmethod
    def adata(self, value: Optional[AnnData]) -> None:
        pass

    @property
    @abstractmethod
    def kernels(self) -> Tuple["KernelExpression", ...]:
        """Underlying base kernels."""

    @property
    @abstractmethod
    def shape(self) -> Tuple[int, int]:
        """``(n_cells, n_cells)``."""

    @property
    @abstractmethod
    def backward(self) -> Optional[bool]:
        """Direction of the process."""

    @abstractmethod
    def __getitem__(self, ix: int) -> "KernelExpression":
        pass

    @abstractmethod
    def __len__(self) -> int:
        pass

    def __add__(self, other: "KernelExpression") -> "KernelExpression":
        return self.__radd__(other)

    def __radd__(self, other: "KernelExpression") -> "KernelExpression":
        def same_level_add(k1: "KernelExpression", k2: "KernelExpression") -> bool:
            if not isinstance(k1, KernelAdd) and not isinstance(k2, KernelMul):
                return False

            for kexpr in k1:
                if not isinstance(kexpr, KernelMul):
                    return False
                if not kexpr._bin_consts:
                    return False

            return True

        if not isinstance(other, KernelExpression):
            return NotImplemented

        s = self * 1.0 if isinstance(self, Kernel) else self
        o = other * 1.0 if isinstance(other, Kernel) else other

        # (c1 * x + c2 * y + ...) + (c3 * z) => (c1 * x + c2 * y + c3 + ... + * z)
        if same_level_add(s, o):
            return KernelAdd(*tuple(s) + (o,))
        if same_level_add(o, s):
            return KernelAdd(*tuple(o) + (s,))

        # add virtual constant
        if not isinstance(s, KernelMul):
            s = s * 1.0
        if not isinstance(o, KernelMul):
            o = o * 1.0

        return KernelAdd(s, o)

    def __mul__(
        self, other: Union[float, int, "KernelExpression"]
    ) -> "KernelExpression":
        return self.__rmul__(other)

    def __rmul__(
        self, other: Union[int, float, "KernelExpression"]
    ) -> "KernelExpression":
        def same_level_mul(k1: "KernelExpression", k2: "KernelExpression") -> bool:
            return (
                isinstance(k1, KernelMul)
                and isinstance(k2, Constant)
                and k1._bin_consts
            )

        if isinstance(other, (int, float, np.integer, np.floating)):
            other = Constant(self.adata, other)

        if not isinstance(other, KernelExpression):
            return NotImplemented

        # fmt: off
        s = self if isinstance(self, (KernelMul, Constant)) else KernelMul(Constant(self.adata, 1.0), self)
        o = other if isinstance(other, (KernelMul, Constant)) else KernelMul(Constant(other.adata, 1.0), other)
        # fmt: on

        # at this point, only KernelMul and Constant is possible (two constants are not possible)
        # (c1 * k) * c2 => (c1 * c2) * k
        if same_level_mul(s, o):
            c, expr = s._bin_const_exprs
            return KernelMul(c * o, expr)
        if same_level_mul(o, s):
            c, expr = o._bin_const_exprs
            return KernelMul(c * s, expr)

        return KernelMul(s, o)

    @d.get_sections(base="write_to_adata", sections=["Parameters"])
    @inject_docs()  # get rid of {{}}
    @d.dedent
    def write_to_adata(self, key: Optional[str] = None) -> None:
        """
        Write the transition matrix and parameters used for computation to the underlying :attr:`adata` object.

        Parameters
        ----------
        key
            Key used when writing transition matrix to :attr:`adata`.
            If `None`, determine the key automatically.

        Returns
        -------
        %(write_to_adata)s
        """
        from cellrank._key import Key

        if self._transition_matrix is None:
            raise ValueError(
                "Compute transition matrix first as `.compute_transition_matrix()`."
            )

        # TODO
        key = Key.uns.kernel(self.backward, key=key)
        # retain the embedding info
        self.adata.uns[f"{key}_params"] = {
            **self.adata.uns.get(f"{key}_params", {}),
            **{"params": self.params},
        }
        self.adata.obsp[key] = self.transition_matrix

    @property
    def transition_matrix(self) -> Union[np.ndarray, csr_matrix]:
        """Row-normalized transition matrix."""
        return self._transition_matrix

    @transition_matrix.setter
    def transition_matrix(self, matrix: Union[np.ndarray, spmatrix]) -> None:
        """
        Set the transition matrix.

        Parameters
        ----------
        matrix
            Transition matrix. If the expression has no parent, the matrix is normalized if needed.

        Returns
        -------
        Nothing, just updates the :attr:`transition_matrix` and optionally normalizes it.
        """
        if matrix.shape != self.shape:
            raise ValueError(
                f"Expected matrix to be of shape `{self.shape}`, found `{matrix.shape}`."
            )

        if issparse(matrix) and not isspmatrix_csr(matrix):
            matrix = csr_matrix(matrix)
        # TODO
        matrix = matrix.astype(np.float64)

        # fmt: off
        should_norm = ~np.isclose(np.asarray(matrix.sum(1)).squeeze(), 1.0, rtol=1e-12).all()
        # fmt: on
        if self._parent is None:
            matrix = _normalize(matrix) if should_norm else matrix
        else:
            matrix = _normalize(matrix) if self._normalize and should_norm else matrix
        # TODO
        if not np.all(np.isfinite(matrix.data if issparse(matrix) else matrix)):
            raise ValueError("Not all values are finite.")

        self._transition_matrix = matrix
        self._params = {}

    @property
    def params(self) -> Dict[str, Any]:
        """Parameters which are used to compute the transition matrix."""
        if len(self.kernels) == 1:
            return self._params
        return {f"{repr(k)}:{i}": k.params for i, k in enumerate(self.kernels)}

    def _format_params(self) -> 3:
        return ", ".join(
            f"{k}={round(v, 3) if isinstance(v, float) else v}"
            for k, v in self.params.items()
        )

    def _reuse_cache(
        self, expected_params: Dict[str, Any], *, time: Optional[Any] = None
    ) -> bool:
        try:
            if expected_params == self._params:
                # TODO
                assert self.transition_matrix is not None
                # assert self.transition_matrix is not None, _ERROR_EMPTY_CACHE_MSG
                # logg.debug(_LOG_USING_CACHE)
                logg.info("    Finish", time=time)
                return True
            return False
        except AssertionError as e:
            raise ValueError(str(e)) from None
        except Exception as e:  # noqa: B902
            logg.error(f"Unable to load the cache, reason: `{e}`")
            return False
        finally:
            self._params = expected_params


class Kernel(KernelExpression, ABC):
    def __init__(self, adata: AnnData, **kwargs: Any):
        super().__init__(**kwargs)
        self._adata = adata
        self._n_obs = adata.n_obs
        self._read_from_adata(**kwargs)

    def _read_from_adata(self, **kwargs: Any) -> None:
        pass

    @property
    def adata(self) -> AnnData:
        return self._adata

    @adata.setter
    def adata(self, adata: Optional[AnnData]) -> None:
        if adata is None:
            self._adata = None
            return
        if not isinstance(adata, AnnData):
            raise TypeError(
                f"Expected `adata` to be of type `AnnData`, found `{type(adata).__name__}`."
            )
        shape = (adata.n_obs, adata.n_obs)
        if self.shape != shape:
            raise ValueError(
                f"Expected the new object to have same shape as the previous `{self.shape}`, found `{shape}`."
            )
        self._adata = adata

    @property
    def kernels(self) -> Tuple["KernelExpression", ...]:
        return (self,)

    @property
    def shape(self) -> Tuple[int, int]:
        return self._n_obs, self._n_obs

    def __getitem__(self, item: int) -> "Kernel":
        return self

    def __len__(self) -> int:
        return 1


class UnidirectionalKernel(UnidirectionalMixin, Kernel, ABC):
    pass


class BidirectionalKernel(BidirectionalMixin, Kernel, ABC):
    pass


class Constant(UnidirectionalKernel):
    def __init__(self, adata: AnnData, value: Union[int, float]):
        super().__init__(adata)
        self.transition_matrix = value

    def compute_transition_matrix(self, value: Union[int, float]) -> "KernelExpression":
        self.transition_matrix = value
        return self

    @property
    def transition_matrix(self) -> Union[int, float]:
        return self._transition_matrix

    @transition_matrix.setter
    def transition_matrix(self, value: Union[int, float]) -> None:
        if not isinstance(value, (int, float, np.integer, np.floating)):
            raise TypeError(
                f"Value must be a `float` or `int`, found `{type(value).__name__}`."
            )
        if value < 0:
            raise ValueError(
                f"Expected the scalar to be non-negative, found `{value}`."
            )

        self._transition_matrix = value
        self._params = {"value": value}

    def copy(self, *, deep: bool = False) -> "KernelExpression":
        return Constant(self.adata, self.transition_matrix)

    def __radd__(
        self, other: Union[int, float, "KernelExpression"]
    ) -> "KernelExpression":
        if isinstance(other, (int, float, np.integer, np.floating)):
            other = Constant(self.adata, other)
        if isinstance(other, Constant):
            if self.shape != other.shape:
                raise ValueError("TODO")
            return Constant(
                self.adata, self.transition_matrix + other.transition_matrix
            )

        return super().__rmul__(other)

    def __rmul__(
        self, other: Union[int, float, "KernelExpression"]
    ) -> "KernelExpression":
        if isinstance(other, (int, float, np.integer, np.floating)):
            other = Constant(self.adata, other)
        if isinstance(other, Constant):
            if self.shape != other.shape:
                raise ValueError("TODO")
            return Constant(
                self.adata, self.transition_matrix * other.transition_matrix
            )

        return super().__rmul__(other)


class NaryKernelExpression(BidirectionalMixin, KernelExpression):
    def __init__(self, *kexprs: KernelExpression):
        super().__init__(parent=None)
        self._validate(kexprs)

    @abstractmethod
    def _combine_transition_matrices(self, t1: Tmat_t, t2: Tmat_t) -> Tmat_t:
        pass

    @property
    @abstractmethod
    def _initial_value(self) -> Union[int, float, Tmat_t]:
        pass

    def compute_transition_matrix(self) -> "KernelExpression":
        for kexpr in self:
            if kexpr.transition_matrix is None:
                if isinstance(kexpr, Kernel):
                    raise RuntimeError(
                        f"Kernel `{kexpr}` is uninitialized. "
                        f"Compute its transition matrix first as `.compute_transition_matrix()`."
                    )
                kexpr.compute_transition_matrix()

        tmat = self._initial_value
        for kexpr in self:
            tmat = self._combine_transition_matrices(tmat, kexpr.transition_matrix)
        self.transition_matrix = tmat

        return self

    def _validate(self, kexprs: Sequence[KernelExpression]) -> None:
        if not len(kexprs):
            raise ValueError("No kernels to combined.")

        shapes = {kexpr.shape for kexpr in kexprs}
        if len(shapes) > 1:
            raise ValueError("TODO")

        directions = {kexpr.backward for kexpr in kexprs}
        if True in directions and False in directions:
            raise ValueError("TODO")

        self._backward = (
            True if True in directions else False if False in directions else None
        )
        self._kexprs = tuple(k if k._parent is None else k.copy() for k in kexprs)
        for kexpr in self:
            kexpr._parent = self

    @property
    def adata(self) -> AnnData:
        return self[0].adata

    @adata.setter
    def adata(self, adata: Optional[AnnData]) -> None:
        # allow re-setting (use for temp. pickling without adata)
        if adata is None or all(kexpr.adata is None for kexpr in self):
            for kexpr in self:
                kexpr.adata = adata
        else:
            self[0].adata = adata

    @property
    def kernels(self) -> Tuple["KernelExpression", ...]:
        kernels = []
        for kexpr in self:
            if isinstance(kexpr, Kernel) and not isinstance(kexpr, Constant):
                kernels.append(kexpr)
            elif isinstance(kexpr, NaryKernelExpression):
                kernels.extend(kexpr.kernels)

        return tuple(kernels)

    def copy(self, *, deep: bool = False) -> "KernelExpression":
        kexprs = tuple(k.copy(deep=deep) for k in self)
        for k in kexprs:
            k._parent = None
        return type(self)(*kexprs)

    @property
    def shape(self) -> Tuple[int, int]:
        return self[0].shape

    def __getitem__(self, ix: int) -> "KernelExpression":
        return self._kexprs[ix]

    def __len__(self) -> int:
        return len(self._kexprs)

    def __invert__(self) -> "KernelExpression":
        if self.backward is None:
            return self.copy()
        kexprs = tuple(
            ~k if isinstance(k, BidirectionalMixin) else k.copy() for k in self
        )
        for k in kexprs:
            k._parent = None
        return type(self)(*kexprs)


class KernelAdd(NaryKernelExpression):
    @property
    def _initial_value(self) -> Union[int, float, Tmat_t]:
        return 0.0

    def _combine_transition_matrices(self, t1: Tmat_t, t2: Tmat_t) -> Tmat_t:
        return t1 + t2


class KernelMul(NaryKernelExpression):
    def compute_transition_matrix(self) -> "KernelExpression":
        self._maybe_recalculate_constants()
        return super().compute_transition_matrix()

    @property
    def _initial_value(self) -> Union[int, float, Tmat_t]:
        return 1.0

    def _combine_transition_matrices(self, t1: Tmat_t, t2: Tmat_t) -> Tmat_t:
        if issparse(t1):
            return t1.multiply(t2)
        if issparse(t2):
            return t2.multiply(t1)
        return t1 * t2

    def _maybe_recalculate_constants(self) -> None:
        constants = self._bin_consts
        if constants:
            total = sum((c.transition_matrix for c in constants), 0.0)
            for c in constants:
                c.transition_matrix = c.transition_matrix / total

            for kexpr in self._kexprs:  # don't normalize  (c * x)
                kexpr._normalize = False

    @property
    def _bin_consts(self) -> List[KernelExpression]:
        if len(self) != 2:
            return []

        return [k for k in self if isinstance(k, Constant)]

    @property
    def _bin_const_exprs(self) -> Tuple[Optional[Constant], Optional[KernelExpression]]:
        if not self._bin_consts:
            return None, None
        k1, k2 = self
        if isinstance(k1, Constant):
            return k1, k2
        return k2, k1
