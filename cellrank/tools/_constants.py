# -*- coding: utf-8 -*-
from enum import Enum
from typing import Union


class PrettyEnum(Enum):
    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self.value)


class Lin(PrettyEnum):
    JOIN = "join"
    REST = "rest"
    NORM = "norm"


class Direction(PrettyEnum):
    FORWARD = "fwd"
    BACKWARD = "bwd"


class RcKey(PrettyEnum):
    FORWARD = "final_states"
    BACKWARD = "root_states"


class LinKey(PrettyEnum):
    FORWARD = "to_final_cells"
    BACKWARD = "from_root_cells"


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
    return f"{k}_cont"


def _dp(k: Union[str, LinKey]) -> str:
    return f"{k}_dp"
