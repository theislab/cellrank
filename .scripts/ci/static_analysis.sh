#!/usr/bin/env bash

pip install black rstcheck
black . --check --diff
rstcheck . --recursive
