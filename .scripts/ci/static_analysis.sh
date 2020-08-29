#!/usr/bin/env bash

set -ev

# make sure that the black version is the same as in .pre-commit-config.yaml, as well as rstcheck
pip install black==20.8b1 git+https://github.com/myint/rstcheck
black . --check --diff --color
rstcheck . --recursive
