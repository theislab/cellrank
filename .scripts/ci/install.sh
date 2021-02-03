#!/usr/bin/env bash

set -ev

python -m pip install --upgrade "pip<21.0.0" --use-deprecated=legacy-resolver
pip install pytest-cov codecov jax jaxlib

if [[ "$OS" == "macos-latest" ]]; then
  pip install -e".[test]" --use-deprecated=legacy-resolver
elif [[ "$OS" == "ubuntu-latest" ]]; then
  if [[ "$USE_SLEPC" == "true" ]]; then
      pip install -e".[krylov,test]"
      python -c "import slepc; import petsc;"
  else
      pip install -e".[test]" --use-deprecated=legacy-resolver
  fi
  pip install rpy2>=3.3.0 --use-deprecated=legacy-resolver
  Rscript --vanilla -e "library('mgcv')"
else
  exit 42
fi
