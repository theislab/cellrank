#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    pip3 install -U pip
    if [[ "$USE_SLEPC" == "true" ]]; then
        pip_cmd=pip3
    fi
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    pip install pytest-cov
    if [[ "$USE_SLEPC" == "true" ]]; then
	echo "Installing SLEPC and PETSc"
        apt-get update
	apt-get install libopenmpi-dev

	pip install mpi4py

	pip install petsc
	pip install petsc4py

	pip install slepc
        pip install slepc4py
    fi
fi
