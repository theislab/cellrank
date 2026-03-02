"""Global CellRank settings."""

import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Final

from rich.console import Console
from rich.logging import RichHandler

__all__ = ["settings"]

_LevelInput = int | str

_LOGGER_NAME: Final[str] = "cellrank"


def _setup_logger() -> logging.Logger:
    """Set up the ``"cellrank"`` logger with a :class:`~rich.logging.RichHandler`."""
    root = logging.getLogger(_LOGGER_NAME)
    if not root.handlers:
        console = Console(stderr=False)
        if console.is_jupyter:
            console.is_jupyter = False
        handler = RichHandler(
            level=logging.DEBUG,
            show_time=False,
            show_path=False,
            rich_tracebacks=True,
            markup=True,
            console=console,
        )
        root.addHandler(handler)
        root.setLevel(logging.INFO)
    root.propagate = False
    return root


class CellRankConfig:
    """Manage user-facing configuration for CellRank."""

    def __init__(self) -> None:
        self._figdir = Path("./figures")

    @property
    def logging_level(self) -> int:
        """Current log level for the ``"cellrank"`` logger."""
        return logging.getLogger(_LOGGER_NAME).level

    @logging_level.setter
    def logging_level(self, level: _LevelInput) -> None:
        """Set the log level for the ``"cellrank"`` logger."""
        logging.getLogger(_LOGGER_NAME).setLevel(level)

    @property
    def figdir(self) -> Path:
        """Default directory for saving figures."""
        return self._figdir

    @figdir.setter
    def figdir(self, value: str | Path) -> None:
        self._figdir = Path(value)

    @contextmanager
    def override_logging_level(self, level: _LevelInput):
        """Temporarily override the logging level.

        Parameters
        ----------
        level
            Level to set for the duration of the context.
        """
        logger = logging.getLogger(_LOGGER_NAME)
        prev = logger.level
        logger.setLevel(level)
        try:
            yield
        finally:
            logger.setLevel(prev)


settings = CellRankConfig()
