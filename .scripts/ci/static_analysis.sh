#!/usr/bin/env bash

set -ev

pip install black==20.8b1 rstcheck
black . --check --diff --color
rstcheck . --recursive
