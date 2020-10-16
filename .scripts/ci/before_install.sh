#!/usr/bin/env bash

set -ev

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    pip3 install -U pip
elif [[ "$TRAVIS_OS_NAME" == "linux" && "$USE_SLEPC" == "true" ]]; then
    echo "Installing APT dependencies"
    sudo apt-get update -y
    sudo apt-get install gcc gfortran libopenmpi-dev libblas-dev liblapack-dev petsc-dev slepc-dev -y
fi

if [[ "$TRAVIS_OS_NAME" == "linux" && "$TRAVIS_BRANCH" == "master" && "$TRAVIS_EVENT_TYPE" == "push" && ! -z "${DEPLOY_TOKEN+x}" ]]; then
    echo "Installing git-lfs"
    sudo apt-get update -y
    sudo apt-get install git-lfs -y
fi
