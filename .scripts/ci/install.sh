#!/usr/bin/env bash

set -e

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    sudo pip3 install -e".[test]"
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    if [[ "$CACHE_NAME" == "krylov" ]]; then
        pip install -e".[krylov,test]"
        python -c "import slepc; import petsc;"
    else
        pip install -e".[test]"
    fi
    pip install pytest-cov
    pip install codecov
fi

python-vendorize
