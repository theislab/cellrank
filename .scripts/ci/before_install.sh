#!/usr/bin/env bash

set -e

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    pip3 install -U pip
elif [[ "$TRAVIS_OS_NAME" == "linux" && "$USE_SLEPC" == "true" ]]; then
    echo "Installing APT dependencies"
    sudo apt-get update -y
    sudo apt-get install gcc gfortran libopenmpi-dev libblas-dev liblapack-dev -y
fi
