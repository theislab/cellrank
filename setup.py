# -*- coding: utf-8 -*-
from setuptools import setup, find_namespace_packages

import os

setup(
    name="cellrank",
    version="1.0.0",
    description="Probabilistic Trajectory Inference using RNA Velocity.",
    url="https://github.com/theislab/cellrank",
    install_requires=list(
        map(str.strip, open(os.path.abspath("requirements.txt"), "r").read().split())
    ),
    zip_safe=False,
    packages=find_namespace_packages(),
    python_required=">=3.6",
)
