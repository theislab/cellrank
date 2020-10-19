#!/usr/bin/env bash

set -ev

sudo apt-get update -y
sudo apt-get install pandoc
pip install -e".[docs,dev]"

pre-commit run --all-files
sphinx-build -j2 -b linkcheck docs/source docs/build/linkcheck
