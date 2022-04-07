#!/usr/bin/env bash

set -euo pipefail

echo "Fetching Yarn GPG key"
# https://github.com/yarnpkg/yarn/issues/7866
curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | sudo apt-key add -

echo "Installing R"
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

echo "Installing APT dependencies"
sudo apt update -y
sudo apt install libopenblas-base r-base r-base-dev r-cran-mgcv -y
