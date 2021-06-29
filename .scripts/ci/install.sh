#!/usr/bin/env bash

set -ev

python -m pip install --upgrade pip
pip install codecov

if [[ "$OS" == "macos-latest" ]]; then
  pip install -e".[test,external]"
elif [[ "$OS" == "ubuntu-latest" ]]; then
  pip install Cython
  if [[ "$USE_SLEPC" == "true" ]]; then
      pip install -e".[krylov,test,external]"
      python -c "import slepc; import petsc;"
  else
      pip install -e".[test,external]"
  fi
  pip install rpy2>=3.3.0
  Rscript --vanilla -e "library('mgcv')"
else
  exit 42
fi
