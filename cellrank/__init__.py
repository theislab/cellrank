# -*- coding: utf-8 -*-
from scanpy import settings

import cellrank.tools as tl
import cellrank.utils as ul
import cellrank.logging
import cellrank.datasets
import cellrank.plotting as pl
from cellrank.tools._read import read
from cellrank.tools._constants import Lin

__author__ = ", ".join(["Marius Lange", "Michal Klein", "Juan Luis Restrepo Lopez"])
__maintainer__ = ", ".join(["Marius Lange", "Michal Klein"])
__email__ = "info@cellrank.org"
__version__ = "1.0.0-rc.2"
