#!/usr/bin/env bash

set -ev

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    sudo pip3 install -e".[test]"
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    if [[ ! -z "${DEPLOY_TOKEN+x}" || "$USE_SLEPC" == "true" ]]; then
        pip install pytest-cov
        pip install codecov
    fi

    if [[ "$TRAVIS_PYTHON_VERSION" == "3.9-dev" ]]; then
        # https://github.com/dhermes/bezier/issues/243
        BEZIER_NO_EXTENSION=true pip install --upgrade bezier --no-binary=bezier
    fi

    if [[ "$USE_SLEPC" == "true" ]]; then
        pip install -e".[krylov,test]"
        python -c "import slepc; import petsc;"
    else
        pip install -e".[test]"
    fi

    pip install rpy2>=3.3.0 jax jaxlib
    Rscript --vanilla -e "library('mgcv')"
fi

python-vendorize
