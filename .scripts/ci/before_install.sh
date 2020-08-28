#!/usr/bin/env bash

set -e

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    pip3 install -U pip
elif [[ "$TRAVIS_OS_NAME" == "linux" && "$USE_SLEPC" == "true" ]]; then
    echo "Installing APT dependencies"
    sudo apt-get update -y
    sudo apt-get install gcc gfortran libopenmpi-dev libblas-dev liblapack-dev -y

    if [[ "$CACHE_NAME" != "krylov" ]]; then
        echo "Installing SLEPc and PETSc Python3 libraries"
        pip_cmd=$(which pip)

        sudo $pip_cmd install mpi4py

        sudo $pip_cmd install petsc
        sudo $pip_cmd install petsc4py

        sudo $pip_cmd install slepc
        sudo $pip_cmd install slepc4py

        python -c "import slepc; import petsc;"
        echo "Successfully installed SLEPc and PETSc"
    fi
fi
