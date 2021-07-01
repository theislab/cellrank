#!/usr/bin/env bash

set -ev

python -m pip install --upgrade pip
# TODO: enable 3.9 CI once new scVelo version is up on PyPI
pip install numpy Cython  # 3.9 external (WOT)
# TODO: add WOT to external once new version is up on PyPI
pip install codecov adjustText git+https://github.com/broadinstitute/wot@master

if [[ "$OS" == "macos-latest" ]]; then
  pip install -e".[test,external]"
elif [[ "$OS" == "ubuntu-latest" ]]; then
  pip install numpy Cython  # 3.9 external (POT)
  if [[ "$USE_SLEPC" == "true" ]]; then
      pip install -e".[krylov,test,external]"
      python -c "import slepc; import petsc;"
  else
      pip install -e".[test,external]"
  fi
  pip install "rpy2>=3.3.0"
  Rscript --vanilla -e "library('mgcv')"
else
  exit 42
fi
