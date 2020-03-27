# -*- coding: utf-8 -*-
from setuptools import setup, find_namespace_packages
from pathlib import Path

import os

try:
    from cellrank import __author__, __email__, __version__
except ImportError:
    __author__ = ""
    __email__ = ""
    __version__ = ""

setup(
    name="cellrank",
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    author=__author__,
    email=__email__,
    maintainer=__author__,
    maintainer_email=__email__,
    version=__version__,
    description="Probabilistic Trajectory Inference using RNA Velocity.",
    long_description=Path("README.rst").read_text("utf-8"),
    url="https://github.com/theislab/cellrank",
    download_url="https://github.com/theislab/cellrank",
    license="BSD",
    install_requires=list(
        map(str.strip, open(os.path.abspath("requirements.txt"), "r").read().split())
    ),
    zip_safe=False,
    packages=find_namespace_packages(),
    python_required=">=3.6",
    platforms=[],
    kewords=[],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Framework :: Jupyter",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux" "Typing :: Typed",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
