from typing import Type, Union, Optional

from abc import abstractmethod
from enum import Enum


class AbstractState(str, Enum):
    @property
    @abstractmethod
    def previous(self) -> Optional["AbstractState"]:
        """TODO."""
        raise NotImplementedError()

    @property
    def format(self) -> str:
        """TODO."""
        return f"Please run `.compute_{self}()` first."

    def raise_exc(self) -> None:
        """TODO."""
        if self.previous is None:
            return
        raise RuntimeError(self.previous.format)

    def __str__(self) -> str:
        return self.value

    def __repr__(self) -> str:
        return str(self)


class StatefulMixin:
    _state_class: Type[AbstractState] = AbstractState

    def __init__(self):
        self._state: Optional[AbstractState] = None

    def _require(self, state: Optional[AbstractState]) -> None:
        if self.state != state:
            # TODO: alt. is to run `state.raise_exc()`
            self.state.raise_exc()

    @property
    def state(self) -> AbstractState:
        return self._state

    @state.setter
    def state(self, value: Union[AbstractState, str]) -> None:
        try:
            self._state = self._state_class(value)
        except ValueError:
            raise ValueError(
                f"Invalid state `{value}`. "
                f"Valid options are `{sorted(self._state_class.__members__.values())}`."
            ) from None
