# TODO: in 2.0, remove .tl
import cellrank.pl
import cellrank.tl
import cellrank.ul
import cellrank.models
import cellrank.kernels
import cellrank.logging
import cellrank.datasets
import cellrank.external
import cellrank.estimators
from cellrank.settings import settings
from cellrank.tl._read import read
from cellrank.tl._lineage import Lin

__author__ = ", ".join(["Marius Lange", "Michal Klein"])
__maintainer__ = ", ".join(["Marius Lange", "Michal Klein"])
__version__ = "1.5.1"
__email__ = "info@cellrank.org"

try:
    from importlib_metadata import version  # Python < 3.8
except ImportError:
    from importlib.metadata import version  # Python = 3.8

from packaging.version import parse

try:
    __full_version__ = parse(version(__name__))
    __full_version__ = (
        f"{__version__}+{__full_version__.local}"
        if __full_version__.local
        else __version__
    )
except ImportError:
    __full_version__ = __version__

del version, parse
