#!/usr/bin/env bash

set -e

pip install black rstcheck
black . --check --diff
rstcheck . --recursive
