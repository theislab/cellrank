#!/usr/bin/env bash

if [[ "$TRAVIS_OS_NAME" == "linux" && "$TRAVIS_PYTHON_VERSION" == "3.8" ]]; then
    pip install codecov
    codecov
fi
