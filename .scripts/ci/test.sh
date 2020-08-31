#!/usr/bin/env bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    python3 -m pytest
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    if [[ ! -z "${DEPLOY_TOKEN+x}" || "$USE_SLEPC" == "true"  ]]; then
        # running regular python -m pytest, even with 1 worker (pytest-parallel or pytest-xdist) causes:
        # MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD with errorcode 50162059.
        python -m pytest --cov-config=.coveragerc --cov=./ --cov-report=xml
    else
        python -m pytest
    fi
fi
