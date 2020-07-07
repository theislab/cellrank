#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    pip3 install -e".[test]"
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    if [[ "$KRYLOV_EXTRA" == "true" ]]; then
        pip install -e.[krylov,test]
    else
        pip install -e.[test]
    fi
fi

python-vendorize
