# -*- coding: utf-8 -*-
import logging
from copy import copy

from scanpy import settings

from cellrank.logging._logging import _RootLogger, _LogFormatter


def _set_log_file(settings):
    file = settings.logfile
    name = settings.logpath
    root = settings._root_logger
    h = logging.StreamHandler(file) if name is None else logging.FileHandler(name)
    h.setFormatter(_LogFormatter())
    h.setLevel(root.level)

    if len(root.handlers) == 1:
        root.removeHandler(root.handlers[0])
    elif len(root.handlers) > 1:
        raise RuntimeError("CellRank's root logger somehow got more than one handler.")

    root.addHandler(h)


settings = copy(settings)
settings._root_logger = _RootLogger(settings.verbosity)
# these 2 lines are necessary to get it working (otherwise no logger is found)
# this is a hacky way of modifying the logging, in the future, use our own
_set_log_file(settings)
settings.verbosity = settings.verbosity
