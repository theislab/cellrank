#!/usr/bin/env bash

set -ev

pip install black rstcheck
black . --check --diff
rstcheck . --recursive
