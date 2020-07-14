#!/usr/bin/env bash

ORIGIN="https://github.com/michalk8/bioconda-recipes"
BIOC_UPSTREAM="https://github.com/bioconda/bioconda-recipes"
PACKAGE="cellrank"

git config --global user.name "TravisCI"

git clone "$ORIGIN"
cd bioconda-recipes

git remote add bioc_upstream "$BIOC_UPSTREAM"
git fetch bioc_upstream
git merge bioc_upstream/master -m "Merge branch 'master' of $UPSTREAM"

# unsetting `PS1` because there are some issues
PS1= ./bootstrap.py /tmp/miniconda
source ~/.config/bioconda/activate

# not including `exclude` causes some issues, just as `no-check-pinnings`
# bioconda seemlingly gets the update on tags pushing (or new PyPI release), no need for `--create-pr`
bioconda-utils autobump recipes/ --packages "$PACKAGE" --no-check-pinnings --exclude ''

git add "recipes/$PACKAGE/meta.yaml"
git commit -m "Update meta.yaml"
git push "https://$DEPLOY_TOKEN@${ORIGIN:8}"
