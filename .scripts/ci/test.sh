#!/usr/bin/env bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    echo "Testing without coverage"
    python3 -m pytest
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    if [[ ! -z "${DEPLOY_TOKEN+x}" ]]; then
        echo "Testing with coverage"
        python -m pytest --cov-config=.coveragerc --cov=./ --cov-report=xml
    else
        echo "Testing without coverage"
        python -m pytest
    fi
fi
