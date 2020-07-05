#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    pip3 install -U pip
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    pip install pytest-cov
    if [[ "$USE_SLEPC" == "true" ]]; then
        echo "Installing dependencies"
        sudo apt-get update -y
        sudo apt-get install gcc gfortran libopenmpi-dev libblas-dev liblapack-dev -y

        pip_cmd=$(which pip)
        echo "Installing SLEPC and PETSc"

        sudo $pip_cmd install mpi4py

        sudo $pip_cmd install petsc
        sudo $pip_cmd install petsc4py

        sudo $pip_cmd install slepc
        sudo $pip_cmd install slepc4py
    fi
fi
