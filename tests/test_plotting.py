# -*- coding: utf-8 -*-
from matplotlib.testing import setup
from pathlib import Path

import cellrank as cr

setup()

HERE: Path = Path(__file__).parent
ROOT = HERE / "_images"
FIGS = HERE / "figures"
DPI = 40

cr.settings.figdir = FIGS
