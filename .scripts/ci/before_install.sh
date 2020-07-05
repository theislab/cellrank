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
	pip install mpi4py --user

	pip install petsc --user
	pip install petsc4py --user

	pip install slepc --user
        pip install slepc4py --user
    fi
fi
