#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    pip3 install -U pip
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    pip install pytest-cov
    if [[ "$USE_SLEPC" == "true" ]]; then
        echo "Installing SLEPC and PETSc"
	      sudo pip install mpi4py

	      sudo pip install petsc
	      sudo pip install petsc4py

	      sudo pip install slepc
        sudo pip install slepc4py
    fi
fi
