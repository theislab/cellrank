#!/usr/bin/env bash

if [[ "$TRAVIS_OS_NAME" == "linux" && ! -z "${DEPLOY_TOKEN+x}" ]]; then
    codecov
fi
