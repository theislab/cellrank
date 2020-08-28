#!/usr/bin/env bash

ORIGIN="https://github.com/michalk8/bioconda-recipes"
DIRNAME="bioconda-recipes"
BIOC_UPSTREAM="https://github.com/bioconda/bioconda-recipes"
PACKAGE="cellrank"


git clone "$ORIGIN" "$DIRNAME"
cd "$DIRNAME"

echo "Fetching $BIOC_UPSTREAM"
git remote add bioc_upstream "$BIOC_UPSTREAM"
git fetch bioc_upstream
git merge bioc_upstream/master -m "Merge branch 'master' of $UPSTREAM"

# unsetting `PS1` because there are some issues
PS1= ./bootstrap.py /tmp/miniconda
source ~/.config/bioconda/activate

# not including `exclude` causes some issues, just as `no-check-pinnings`
# bioconda seemlingly gets the update on tags pushing (or new PyPI release), no need for `--create-pr`
echo "Updating meta.yaml"
bioconda-utils autobump recipes/ --packages "$PACKAGE" --no-check-pinnings --exclude ''

git add -f "recipes/$PACKAGE/meta.yaml"
git commit -m "[ci skip] Update meta.yaml: $TRAVIS_BUILD_ID"

if [[ $? -eq 0 ]]; then
    echo "Pushing changes"
    git push "https://$DEPLOY_TOKEN@${ORIGIN:8}" --quiet
else
    echo "Nothing to commit, this shouldn't have happened"
fi

cd ..
rm -rf "$DIRNAME"
