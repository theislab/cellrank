#!/usr/bin/env bash

set -ev

pip install pre-commit
# this doesn't show which things would be blackened, but it should be fine
pre-commit run --all-files
