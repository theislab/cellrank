#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    pip3 install -e".[test]"
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    pip install -e.[test]
fi

python-vendorize
