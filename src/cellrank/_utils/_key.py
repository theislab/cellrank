from typing import Any, Callable, Optional

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
    def backward(cls, bwd: Optional[bool]) -> str:
        return "bwd" if bwd else "fwd"

    @classmethod
    def where(cls, bwd: Optional[bool]) -> str:
        return "from" if bwd else "to"

    @classmethod
    def initial(cls, bwd: Optional[bool]) -> str:
        return "initial" if bwd else "terminal"

    @classmethod
    def cytotrace(cls, key: str) -> str:
        return f"ct_{key}"

    class obs:
        @classmethod
        def probs(cls, key: str) -> str:
            return f"{key}_probs"

        @classmethod
        def macrostates(cls, bwd: Optional[bool]) -> str:
            return f"macrostates_{Key.backward(bwd)}"

        @classmethod
        def term_states(cls, estim_bwd: Optional[bool], *, bwd: bool = False) -> str:
            states = "init_states" if bwd else "term_states"
            return f"{states}_{Key.backward(estim_bwd)}"

        @classmethod
        def priming_degree(cls, bwd: Optional[bool]) -> str:
            return f"priming_degree_{Key.backward(bwd)}"

    class obsm:
        @classmethod
        def memberships(cls, key: str) -> str:
            return f"{key}_memberships"

        @classmethod
        def schur_vectors(cls, bwd: Optional[bool]) -> str:
            return f"schur_vectors_{Key.backward(bwd)}"

        @classmethod
        def macrostates(cls, bwd: Optional[bool]) -> str:
            return f"macrostates_{Key.backward(bwd)}"

        @classmethod
        def fate_probs(cls, bwd: Optional[bool]) -> str:
            return f"lineages_{Key.backward(bwd)}"

        @classmethod
        def abs_times(cls, bwd: Optional[bool]) -> str:
            return f"absorption_times_{Key.backward(bwd)}"

    class varm:
        @classmethod
        def lineage_drivers(cls, bwd: Optional[bool]):
            return Key.initial(bwd) + "_lineage_drivers"

    class uns:
        @classmethod
        def kernel(cls, bwd: Optional[bool], key: Optional[str] = None) -> str:
            return key if key is not None else f"T_{Key.backward(bwd)}"

        @classmethod
        def estimator(cls, bwd: Optional[bool], key: Optional[str] = None) -> str:
            return key if key is not None else f"{Key.backward(bwd)}_estimator"

        @classmethod
        def names(cls, key: str) -> str:
            return f"{key}_names"

        @classmethod
        def colors(cls, key: str) -> str:
            return f"{key}_colors"

        @classmethod
        def eigen(cls, bwd: Optional[bool]) -> str:
            return f"eigendecomposition_{Key.backward(bwd)}"

        @classmethod
        def schur_matrix(cls, bwd: Optional[bool]) -> str:
            return f"schur_matrix_{Key.backward(bwd)}"

        @classmethod
        def coarse(cls, bwd: Optional[bool]) -> str:
            return f"coarse_{Key.backward(bwd)}"
