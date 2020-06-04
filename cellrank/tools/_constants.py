# -*- coding: utf-8 -*-
"""Module containing CellRank contants."""

from enum import Enum
from typing import Union


class PrettyEnum(Enum):
    """Enum wit a pretty __str__ and __repr__."""

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.value)


class Lin(PrettyEnum):
    """Lineage aggregation type."""

    REST = "rest"
    OTHERS = "others"


class Direction(PrettyEnum):
    """Direction of the process."""

    FORWARD = "fwd"
    BACKWARD = "bwd"


class StateKey(PrettyEnum):
    """State key in `adata.obs`."""

    FORWARD = "final_states"
    BACKWARD = "root_states"


# StateKey and LinKey must have the same suffix `_..._states` because of model.prepare
class LinKey(PrettyEnum):
    """Lineage key in `adata.obsm`."""

    FORWARD = "to_final_states"
    BACKWARD = "from_root_states"


class MetaKey(PrettyEnum):
    """Metastable state key in `adata.obs`."""

    FORWARD = "metastable_states_fwd"
    BACKWARD = "metastable_states_bwd"


class Prefix(PrettyEnum):
    """Direction prefix."""

    FORWARD = "to"
    BACKWARD = "from"


def _transition(d: Union[str, Direction]) -> str:
    return f"T_{d}"


def _lin_names(k: Union[str, LinKey]) -> str:
    return f"{k}_names"


def _colors(k: Union[str, LinKey, StateKey]) -> str:
    return f"{k}_colors"


def _probs(k: Union[str, StateKey]) -> str:
    return f"{k}_probs"


def _dp(k: Union[str, LinKey]) -> str:
    return f"{k}_dp"
