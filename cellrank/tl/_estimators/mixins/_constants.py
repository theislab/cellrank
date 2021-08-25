from typing import Any, Callable


class cprop:
    def __init__(self, f: Callable[..., str]):
        self.f = f

    def __get__(self, obj: Any, owner: Any) -> str:
        return self.f(owner)


class Key:
    @classmethod
    def backward(cls, bwd: bool) -> str:
        return "bwd" if bwd else "fwd"

    @classmethod
    def where(cls, bwd: bool) -> str:
        return "from" if bwd else "to"

    class obs:
        @classmethod
        def probs(self, key: str) -> str:
            return f"{key}_probs"

        @classmethod
        def term_states(cls, bwd: bool) -> str:
            return ("initial" if bwd else "terminal") + "_states"

        @classmethod
        def priming_degree(cls, bwd: bool) -> str:
            return f"priming_degree_{Key.backward(bwd)}"

    class obsm:
        @classmethod
        def macrostates(cls, bwd: bool) -> str:
            return f"macrostates_{Key.backward(bwd)}"

        @classmethod
        def abs_probs(cls, bwd: bool) -> str:
            return ("from" if bwd else "to") + "_" + Key.obs.term_states(bwd)

        @classmethod
        def abs_times(cls, bwd: bool) -> str:
            return f"absorption_times_{Key.backward(bwd)}"

    class varm:
        @classmethod
        def lineage_drivers(cls, bwd: bool):
            return ("initial" if bwd else "terminal") + "_lineage_drivers"

    class uns:
        @classmethod
        def names(cls, key: str) -> str:
            return f"{key}_names"

        @classmethod
        def colors(cls, key: str) -> str:
            return f"{key}_colors"

        @classmethod
        def eigen(cls, bwd: bool) -> str:
            return f"eigendecomposition_{Key.backward(bwd)}"
