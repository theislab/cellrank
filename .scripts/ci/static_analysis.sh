#!/usr/bin/env bash

set -ev

pip install -e".[docs,dev]"
# this doesn't show which things would be blackened, but it should be fine
pre-commit run --all-files
sphinx-build -j2 -b linkcheck docs/source docs/build/linkcheck
