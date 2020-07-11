#!/usr/bin/env bash

UPSTREAM="https://github.com/bioconda/bioconda-recipes"

git clone https://github.com/michalk8/bioconda-recipes
cd bioconda-recipes

git remote add upstream "$UPSTREAM"
git fetch upstream
git merge upstream/master -m "Merge branch 'master' of $UPSTREAM"

# unsetting `PS1` because there are some issues
PS1= ./bootstrap.py /tmp/miniconda
source ~/.config/bioconda/activate

# not including `exclude` causes some issue, just as `no-check-pinnings`
GITHUB_TOKEN="$DEPLOY_TOKEN" bioconda-utils autobump recipes/ --packages cellrank --no-check-pinnings --exclude '' --create-pr
