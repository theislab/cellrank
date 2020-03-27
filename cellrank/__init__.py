# -*- coding: utf-8 -*-
import cellrank.tools as tl
import cellrank.plotting as pl
import cellrank.utils as ul
import cellrank.datasets

from scanpy import settings


__author__ = "".join([])
__email__ = "".join([])

try:
    from setuptools_scm import get_version

    __version__ = get_version(root="..", relative_to=__file__)
except (LookupError, ImportError):
    from packaging import version

    try:
        from importlib_metadata import version as get_version  # Python < 3.8
    except ImportError:
        from importlib.metadata import version as get_version  # Python = 3.8

    __version__ = version.parse(get_version(__name__))
    del version

del get_version
