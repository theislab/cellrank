# -*- coding: utf-8 -*-
from enum import Enum
from typing import Union


class PrettyEnum(Enum):
    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.value)


class Lin(PrettyEnum):
    REST = "rest"
    OTHERS = "others"


class Direction(PrettyEnum):
    FORWARD = "fwd"
    BACKWARD = "bwd"


class RcKey(PrettyEnum):
    FORWARD = "final_states"
    BACKWARD = "root_states"


# RcKey and LinKey must have the same suffix `_..._states` because of model.prepare
class LinKey(PrettyEnum):
    FORWARD = "to_final_states"
    BACKWARD = "from_root_states"


class MetaKey(PrettyEnum):
    FORWARD = "metastable_states_fwd"
    BACKWARD = "metastable_states_bwd"


class Prefix(PrettyEnum):
    FORWARD = "to"
    BACKWARD = "from"


def _transition(d: Union[str, Direction]) -> str:
    return f"T_{d}"


def _lin_names(k: Union[str, LinKey]) -> str:
    return f"{k}_names"


def _colors(k: Union[str, LinKey, RcKey]) -> str:
    return f"{k}_colors"


def _probs(k: Union[str, RcKey]) -> str:
    return f"{k}_probs"


def _dp(k: Union[str, LinKey]) -> str:
    return f"{k}_dp"
