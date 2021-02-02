#!/usr/bin/env bash

set -ev

python -m pip install --upgrade pip
pip install pytest-cov codecov jax jaxlib

if [[ "$OS" == "macos-latest" ]]; then
  pip install -e".[test]"
elif [[ "$OS" == "ubuntu-latest" ]]; then
  if [[ "$USE_SLEPC" == "true" ]]; then
      pip install -e".[krylov,test]"
      python -c "import slepc; import petsc;"
  else
      pip install -e".[test]"
  fi
  pip install rpy2>=3.3.0
  Rscript --vanilla -e "library('mgcv')"
else
  exit 42
fi
