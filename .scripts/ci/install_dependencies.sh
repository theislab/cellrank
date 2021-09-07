#!/usr/bin/env bash

set -ev

if [[ "$RUNNER_OS" == "Linux" ]]; then
    echo "Installing APT dependencies"
    # https://github.com/yarnpkg/yarn/issues/7866
    curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | sudo apt-key add -
    # R-related
    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

    sudo apt update -y
    sudo apt install libopenblas-base r-base r-base-dev r-cran-mgcv -y

    if [[ "$USE_SLEPC" == "true" ]]; then
        echo "Installing PETSc/SLEPc dependencies"
        sudo apt install gcc gfortran libopenmpi-dev libblas-dev liblapack-dev petsc-dev slepc-dev -y
    fi
elif [[ "$RUNNER_OS" != "macOS" ]]; then
    echo "Invalid OS for dependencie installations: $OS"
    exit 42
fi
