#!/usr/bin/env bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    python3 -m pytest
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    pip install pytest-cov
    python -m pytest --cov-config=.coveragerc --cov=./ --cov-report=xml
fi
