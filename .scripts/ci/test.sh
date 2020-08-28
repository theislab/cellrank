#!/usr/bin/env bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    python3 -m pytest
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    python -m pytest --cov-config=.coveragerc --cov=./ --cov-report=xml
else
    echo "Windows build is not supported"
    exit 1
fi
