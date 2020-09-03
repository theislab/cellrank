# -*- coding: utf-8 -*-
from scanpy import settings

import cellrank.pl
import cellrank.tl
import cellrank.ul
import cellrank.logging
import cellrank.datasets
from cellrank.tl._read import read
from cellrank.tl._constants import Lin

__author__ = ", ".join(["Marius Lange", "Michal Klein", "Juan Luis Restrepo Lopez"])
__maintainer__ = ", ".join(["Marius Lange", "Michal Klein"])
__version__ = "1.0.0-rc.6"
__email__ = "info@cellrank.org"
