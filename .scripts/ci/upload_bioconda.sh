#!/usr/bin/env bash

git clone https://github.com/michalk8/bioconda-recipes
cd bioconda-recipes

git remote add upstream "$UPSTREAM"
git fetch upstream
git merge upstream/master -m "Merge branch 'master' of $UPSTREAM"

# TODO update version + hash

git add recipes/cellrank/meta.yaml
git commit -m "Update meta.yaml"
git push

PS1= ./bootstrap.py /tmp/miniconda
source ~/.config/bioconda/activate

bioconda-utils autobump recipes/ --packages cellrank --no-check-pinnings --exclude ''
