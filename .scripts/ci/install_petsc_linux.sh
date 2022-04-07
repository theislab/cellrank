#!/usr/bin/env bash

set -euo pipefail

function install_petsc {
    echo "Installing PETSc; version: '$PC_VERSION', dir: '$PETSC_DIR' with arch: '$PETSC_ARCH'"
    curl -O "https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-$PC_VERSION.tar.gz"
    tar -xzf "petsc-lite-$PC_VERSION.tar.gz"
    pushd "petsc-$PC_VERSION"
    ./configure --with-cc=mpicc --with-fc=0 --with-cxx=mpicxx --with-debugging=0 --with-mpi=1
    make all
    make check
    popd
}

function install_slepc {
    echo "Installing SLEPc; version: '$SC_VERSION', dir: '$SLEPC_DIR'"
    curl -O "https://slepc.upv.es/download/distrib/slepc-$SC_VERSION.tar.gz"
    tar -xzf "slepc-$SC_VERSION.tar.gz"
    pushd "slepc-$SC_VERSION"
    ./configure
    make all
    make check
    popd
}

echo "Installing PETSc/SLEPc dependencies"
sudo apt-get update -y
sudo apt-get install gcc gfortran libopenmpi-dev libblas-dev liblapack-dev -y

pushd "$HOME"
install_petsc
install_slepc
popd

echo "Installing numpy"
python -m pip install --upgrade pip
pip install numpy

echo "Symlinking numpy"
NUMPY_INCLUDE="$(python -c 'import numpy; print(numpy.get_include())')"
ln -sfv "$NUMPY_INCLUDE/numpy" "$PETSC_DIR/$PETSC_ARCH/include"
