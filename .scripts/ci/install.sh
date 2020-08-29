#!/usr/bin/env bash

set -ev

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    sudo pip3 install -e".[test]"
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    if [[ ! -z "${DEPLOY_TOKEN+x}" || "$USE_SLEPC" == "true" ]]; then
        pip install pytest-cov
        pip install codecov
    fi
    if [[ "$USE_SLEPC" == "true" ]]; then
        pip install -e".[krylov,test]"
        python -c "import slepc; import petsc;"
    else
        pip install -e".[test]"
    fi
fi

python-vendorize
