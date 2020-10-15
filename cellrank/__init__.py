# -*- coding: utf-8 -*-
import cellrank.pl
import cellrank.tl
import cellrank.ul
import cellrank.logging
import cellrank.datasets
from cellrank.settings import settings as settings
from cellrank.tl._read import read
from cellrank.tl._constants import Lin

__author__ = ", ".join(["Marius Lange", "Michal Klein", "Juan Luis Restrepo Lopez"])
__maintainer__ = ", ".join(["Marius Lange", "Michal Klein"])
__version__ = "1.0.0-rc.12"
__email__ = "info@cellrank.org"

try:
    from importlib_metadata import version  # Python < 3.8
except ImportError:
    from importlib.metadata import version  # Python = 3.8
__full_version__ = version(__name__)

del version
