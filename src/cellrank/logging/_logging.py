import datetime
import functools
import logging
from typing import Iterable, Optional

__all__ = [
    "print_versions",
    "print_version_and_date",
    "hint",
    "info",
    "warning",
    "error",
    "debug",
]

HINT = (logging.INFO + logging.DEBUG) // 2
logging.addLevelName(HINT, "HINT")


class _RootLogger(logging.RootLogger):
    def __init__(self, level):
        super().__init__(level)
        self.propagate = False
        _RootLogger.manager = logging.Manager(self)

    def log(
        self,
        level: int,
        msg: str,
        *,
        extra: Optional[dict] = None,
        time: datetime.datetime = None,
        deep: Optional[str] = None,
    ) -> datetime.datetime:
        from cellrank import settings

        # this will correctly initialize the handles if doing
        # just from cellrank import logging
        settings.verbosity = settings.verbosity

        now = datetime.datetime.now(datetime.timezone.utc)
        time_passed: datetime.timedelta = None if time is None else now - time
        extra = {
            **(extra or {}),
            "deep": deep if settings.verbosity.level < level else None,
            "time_passed": time_passed,
        }
        _ = super().log(level, msg, extra=extra)
        return now

    def critical(self, msg: str, *, time=None, deep=None, extra=None) -> datetime.datetime:
        return self.log(logging.CRITICAL, msg, time=time, deep=deep, extra=extra)

    def error(self, msg: str, *, time=None, deep=None, extra=None) -> datetime.datetime:
        return self.log(logging.ERROR, msg, time=time, deep=deep, extra=extra)

    def warning(self, msg: str, *, time=None, deep=None, extra=None) -> datetime.datetime:
        return self.log(logging.WARNING, msg, time=time, deep=deep, extra=extra)

    def info(self, msg: str, *, time=None, deep=None, extra=None) -> datetime.datetime:
        return self.log(logging.INFO, msg, time=time, deep=deep, extra=extra)

    def hint(self, msg: str, *, time=None, deep=None, extra=None) -> datetime.datetime:
        return self.log(HINT, msg, time=time, deep=deep, extra=extra)

    def debug(self, msg: str, *, time=None, deep=None, extra=None) -> datetime.datetime:
        return self.log(logging.DEBUG, msg, time=time, deep=deep, extra=extra)


class _LogFormatter(logging.Formatter):
    def __init__(self, fmt="{levelname}: {message}", datefmt="%Y-%m-%d %H:%M", style="{"):
        super().__init__(fmt, datefmt, style)

    def format(self, record: logging.LogRecord):
        format_orig = self._style._fmt
        if record.levelno == logging.INFO:
            self._style._fmt = "{message}"
        elif record.levelno == HINT:
            self._style._fmt = "--> {message}"
        elif record.levelno == logging.DEBUG:
            self._style._fmt = "DEBUG: {message}"
        if record.time_passed:
            # strip microseconds
            if record.time_passed.microseconds:
                record.time_passed = datetime.timedelta(seconds=int(record.time_passed.total_seconds()))
            if "{time_passed}" in record.msg:
                record.msg = record.msg.replace("{time_passed}", str(record.time_passed))
            else:
                self._style._fmt += " ({time_passed})"
        if record.deep:
            record.msg = f"{record.msg}: {record.deep}"
        result = logging.Formatter.format(self, record)
        self._style._fmt = format_orig
        return result


_DEPENDENCIES_NUMERICS = [
    "scanpy",
    "anndata",
    "numpy",
    "numba",
    "scipy",
    "pandas",
    "pygpcca",
    ("sklearn", "scikit-learn"),
    "statsmodels",
    ("igraph", "python-igraph"),
    "scvelo",
    "pygam",
]


_DEPENDENCIES_PLOTTING = ["matplotlib", "seaborn"]


def _versions_dependencies(dependencies: Iterable[str]):
    # this is not the same as the requirements!
    for mod in dependencies:
        mod_name, dist_name = mod if isinstance(mod, tuple) else (mod, mod)
        try:
            imp = __import__(mod_name)
            yield dist_name, imp.__version__
        except (ImportError, AttributeError):
            pass


def print_versions():
    """Print package versions that might influence the numerical and plotting results."""
    from cellrank import settings

    modules = ["cellrank"] + _DEPENDENCIES_NUMERICS + _DEPENDENCIES_PLOTTING
    print(
        " ".join(f"{mod}=={ver}" for mod, ver in _versions_dependencies(modules)),
        file=settings.logfile,
    )


def print_version_and_date():
    """
    Print version and date.

    Useful for starting a notebook so you see when you started working.
    """
    from cellrank import __version__, settings

    print(
        f"Running CellRank {__version__}, on {datetime.datetime.now():%Y-%m-%d %H:%M}.",
        file=settings.logfile,
    )


def _copy_docs_and_signature(fn):
    return functools.partial(functools.update_wrapper, wrapped=fn, assigned=["__doc__", "__annotations__"])


def error(
    msg: str,
    *,
    time: datetime.datetime = None,
    deep: Optional[str] = None,
    extra: Optional[dict] = None,
) -> datetime.datetime:
    """
    Log message with specific level and return current time.

    Parameters
    ----------
    msg
        Message to display.
    time
        A time in the past. If this is passed, the time difference from then
        to now is appended to `msg` as ` (HH:MM:SS)`.
        If `msg` contains `{time_passed}`, the time difference is instead
        inserted at that position.
    deep
        If the current verbosity is higher than the log function’s level,
        this gets displayed as well
    extra
        Additional values you can specify in `msg` like `{time_passed}`.

    Returns
    -------
    The current time.
    """
    from cellrank import settings

    return settings._root_logger.error(msg, time=time, deep=deep, extra=extra)


@_copy_docs_and_signature(error)
def warning(msg: str, *, time=None, deep=None, extra=None) -> datetime.datetime:
    from cellrank import settings

    return settings._root_logger.warning(msg, time=time, deep=deep, extra=extra)


@_copy_docs_and_signature(error)
def info(msg: str, *, time=None, deep=None, extra=None) -> datetime.datetime:
    from cellrank import settings

    return settings._root_logger.info(msg, time=time, deep=deep, extra=extra)


@_copy_docs_and_signature(error)
def hint(msg: str, *, time=None, deep=None, extra=None) -> datetime.datetime:
    from cellrank import settings

    return settings._root_logger.hint(msg, time=time, deep=deep, extra=extra)


@_copy_docs_and_signature(error)
def debug(msg: str, *, time=None, deep=None, extra=None) -> datetime.datetime:
    from cellrank import settings

    return settings._root_logger.debug(msg, time=time, deep=deep, extra=extra)
