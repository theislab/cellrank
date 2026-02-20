from collections.abc import Callable
from typing import Any

__all__ = ["Key"]


class cprop:
    """Class property."""

    def __init__(self, f: Callable[..., str]):
        self.f = f

    def __get__(self, obj: Any, owner: Any) -> str:
        return self.f(owner)


class Key:
    """Class which manages keys in :class:`anndata.AnnData`."""

    @classmethod
    def backward(cls, bwd: bool | None) -> str:
        return "bwd" if bwd else "fwd"

    @classmethod
    def where(cls, bwd: bool | None) -> str:
        return "from" if bwd else "to"

    @classmethod
    def initial(cls, bwd: bool | None) -> str:
        return "initial" if bwd else "terminal"

    @classmethod
    def cytotrace(cls, key: str) -> str:
        return f"ct_{key}"

    class obs:
        @classmethod
        def probs(cls, key: str) -> str:
            return f"{key}_probs"

        @classmethod
        def macrostates(cls, bwd: bool | None) -> str:
            return f"macrostates_{Key.backward(bwd)}"

        @classmethod
        def term_states(cls, estim_bwd: bool | None, *, bwd: bool = False) -> str:
            states = "init_states" if bwd else "term_states"
            return f"{states}_{Key.backward(estim_bwd)}"

        @classmethod
        def priming_degree(cls, bwd: bool | None) -> str:
            return f"priming_degree_{Key.backward(bwd)}"

    class obsm:
        @classmethod
        def memberships(cls, key: str) -> str:
            return f"{key}_memberships"

        @classmethod
        def schur_vectors(cls, bwd: bool | None) -> str:
            return f"schur_vectors_{Key.backward(bwd)}"

        @classmethod
        def macrostates(cls, bwd: bool | None) -> str:
            return f"macrostates_{Key.backward(bwd)}"

        @classmethod
        def fate_probs(cls, bwd: bool | None) -> str:
            return f"lineages_{Key.backward(bwd)}"

        @classmethod
        def abs_times(cls, bwd: bool | None) -> str:
            return f"absorption_times_{Key.backward(bwd)}"

    class varm:
        @classmethod
        def lineage_drivers(cls, bwd: bool | None):
            return Key.initial(bwd) + "_lineage_drivers"

    class uns:
        @classmethod
        def kernel(cls, bwd: bool | None, key: str | None = None) -> str:
            return key if key is not None else f"T_{Key.backward(bwd)}"

        @classmethod
        def estimator(cls, bwd: bool | None, key: str | None = None) -> str:
            return key if key is not None else f"{Key.backward(bwd)}_estimator"

        @classmethod
        def names(cls, key: str) -> str:
            return f"{key}_names"

        @classmethod
        def colors(cls, key: str) -> str:
            return f"{key}_colors"

        @classmethod
        def eigen(cls, bwd: bool | None) -> str:
            return f"eigendecomposition_{Key.backward(bwd)}"

        @classmethod
        def schur_matrix(cls, bwd: bool | None) -> str:
            return f"schur_matrix_{Key.backward(bwd)}"

        @classmethod
        def coarse(cls, bwd: bool | None) -> str:
            return f"coarse_{Key.backward(bwd)}"
