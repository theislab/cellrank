#!/usr/bin/env bash

ORIGIN="https://github.com/theislab/cellrank_notebooks"
DIRNAME="cellrank_notebooks"

if [[ -n "${DEPLOY_TOKEN+x}" ]]; then
  git clone "$ORIGIN" "$DIRNAME"
  cd cellrank_notebooks

  cp ../README.rst .

  git add "README.rst"
  git commit -m "Update README.rst"  # don't use [ci skip], we want to regenerate the notebooks
  git push "https://$DEPLOY_TOKEN@${ORIGIN:8}"

  cd ..
  rm -rf "$DIRNAME"

fi
