#!/usr/bin/env bash

set -ev

if [[ "$OS" == "ubuntu-latest" ]]; then
    echo "Installing APT dependencies"
    # R-related
    sudo apt-key adv --recv-keys --keyserver "http://keyserver.ubuntu.com:11371/pks/lookup?op=get&search=0xE298A3A825C0D65DFD57CBB651716619E084DAB9"
    sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/"

    sudo apt-get update -y
    sudo apt-get install libopenblas-base r-base r-base-dev r-cran-mgcv -y
  if [[ "$USE_SLEPC" == "true" ]]; then
    echo "Installing PETSc SLEPc dependencies"
    sudo apt-get install gcc gfortran libopenmpi-dev libblas-dev liblapack-dev petsc-dev slepc-dev -y
  fi
elif [[ "$OS" != "macos-latest" ]]; then
  exit 42
fi
