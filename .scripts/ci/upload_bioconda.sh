#!/usr/bin/env bash

ORIGIN="https://github.com/michalk8/bioconda-recipes"
UPSTREAM="https://github.com/bioconda/bioconda-recipes"

git clone "$ORIGIN"
cd bioconda-recipes

git remote add upstream "$UPSTREAM"
git fetch upstream
git merge upstream/master -m "Merge branch 'master' of $UPSTREAM"

# unsetting `PS1` because there are some issues
PS1= ./bootstrap.py /tmp/miniconda
source ~/.config/bioconda/activate

# not including `exclude` causes some issue, just as `no-check-pinnings`
GITHUB_TOKEN="$DEPLOY_TOKEN" bioconda-utils autobump recipes/ --packages cellrank --no-check-pinnings --exclude '' --create-pr

git add recipes/cellrank/meta.yaml
git push "https://$DEPLOY_TOKEN@${ORIGIN:8}" -m "Update CellRank meta.yaml"
