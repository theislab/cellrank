![PyPI](https://img.shields.io/pypi/v/cellrank)
[![Build Status](https://travis-ci.com/theislab/cellrank.svg?token=UqaQZpSKCK4qZNfu1sqq&branch=master)](https://travis-ci.com/theislab/cellrank)
![GitHub](https://img.shields.io/github/license/theislab/cellrank)

# CellRank
Probabilistic Trajectory Inference using RNA Velocity

## Notebooks
Tutorial notebooks are in a separate repository [here](https://github.com/theislab/cellrank_notebooks).


## Installation
```bash
git clone https://github.com/theislab/cellrank.git
cd cellrank && pip install .
```
or
```bash
pip install git+https://github.com/theislab/cellrank
```

## Documentation
We already have a sphinx build documentation but we can't host it yet becasue this repo is still private. So for now, if you would like to take a look at the documentation, install [sphinx](https://www.sphinx-doc.org/en/master/usage/installation.html), navigate into the `docs` folder and run `make html` in your terminal, which will generate a set of `html` pages in the `build` folder. Just click on `index.html` to see the documentation. For more info, see the [sphinx quickstart guide](https://www.sphinx-doc.org/en/master/usage/quickstart.html)

## Contributions
We welcome your contributions to CellRank! Please make sure to have a look at the guidelines below:
- git commits: we follow [these](https://chris.beams.io/posts/git-commit/) rules
- docstrings: we follow the numpydoc style in general and (scanpy's style)[https://github.com/theislab/scanpy/blob/master/CONTRIBUTING.md] in particular

## Gene Trend Plotting
We observed that using R libraries for gene trend plotting gives us better results, therefore, the default is to use these librarires via `rpy2`. If you don't have `rpy2` installed, we use `pygam` instead. If you use anaconda or miniconda, consider installing `rpy2` via `conda install -c r rpy2`, see [here](https://anaconda.org/r/rpy2).
