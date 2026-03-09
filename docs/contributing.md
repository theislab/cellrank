# Contributing guide

This document summarizes the most important information for getting you started
on contributing to CellRank.
We assume you are already familiar with git and with making pull requests on GitHub.
For more extensive tutorials covering the absolute basics, refer to the
[pyOpenSci tutorials], the [Scientific Python tutorials], or the [scanpy developer guide].

[pyOpenSci tutorials]: https://www.pyopensci.org/learn.html
[Scientific Python tutorials]: https://learn.scientific-python.org/development/tutorials/
[scanpy developer guide]: https://scanpy.readthedocs.io/en/latest/dev/index.html

:::{tip} The *hatch* project manager
We recommend familiarizing yourself with [hatch].
Hatch is a Python project manager that

- manages virtual environments, separately for development, testing, and building the documentation.
  Separating environments avoids dependency conflicts.
- allows you to run tests locally in different environments (e.g., different Python versions).
- allows you to run tasks defined in `pyproject.toml`, e.g., to build documentation.

While the project is set up with `hatch` in mind,
you can also use other tools such as `uv` or `pip`.
:::

[hatch]: https://hatch.pypa.io/latest/

## Installing dev dependencies

In addition to the packages needed to _use_ CellRank,
you need additional packages to [run tests](#writing-tests) and [build the documentation](#building-the-docs-locally).

:::::{tab-set}
::::{tab-item} Hatch
:sync: hatch

On the command line, you typically interact with hatch through its CLI.
Running one of the following commands will automatically resolve the environments
for testing and building the documentation in the background:

```bash
hatch test  # defined in [tool.hatch.envs.hatch-test] in pyproject.toml
hatch run docs:build  # defined in [tool.hatch.envs.docs]
```

When using an IDE such as VS Code,
you'll have to point the editor at the virtual environment paths manually.
The environment you typically want to use is the `hatch-test`
environment with the latest Python version.

To get a list of all environments for your project, run

```bash
hatch env show -i
```

Select the environment name you want from the `Envs` column
(e.g., `hatch-test.py3.12-stable`) and create it:

```bash
hatch env create hatch-test.py3.12-stable
```

Then, obtain the path to the environment using

```bash
hatch env find hatch-test.py3.12-stable
```

In VS Code, open the command palette (`Ctrl+Shift+P` / `Cmd+Shift+P`) and search
for `Python: Select Interpreter`. Choose `Enter Interpreter Path` and paste the
path from above.
::::

::::{tab-item} uv
:sync: uv

A popular choice for managing virtual environments is [uv].
The main disadvantage compared to hatch is that it supports only a single
environment per project at a time, which requires mixing dependencies for
running tests and building docs.

To initialize a virtual environment in the `.venv` directory, simply run

```bash
uv sync --all-extras
```

The `.venv` directory is typically automatically discovered by IDEs such as VS Code.
::::

::::{tab-item} Pip
:sync: pip

You can also manage environments manually using `pip`:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,test,doc]"
```

The `.venv` directory is typically automatically discovered by IDEs such as VS Code.
::::
:::::

[uv]: https://docs.astral.sh/uv/

## Code-style

This package uses [pre-commit] to enforce consistent code-styles.
On every commit, pre-commit checks will either automatically fix issues
with the code, or raise an error message.

To enable pre-commit locally, run

```bash
pre-commit install
```

in the root of the repository.
Pre-commit will automatically download all dependencies when it is run for the first time.

Alternatively, you can rely on the [pre-commit.ci] service enabled on GitHub.
If you didn't run `pre-commit` before pushing changes to GitHub, it will
automatically commit fixes to your pull request, or show an error message.
If pre-commit.ci added a commit on a branch you still have been working on locally, use

```bash
git pull --rebase
```

to integrate the changes into yours.
While the [pre-commit.ci] is useful, we strongly encourage installing and running pre-commit
locally first to understand its usage.

We use [Ruff] for linting and formatting Python code,
and [Biome] for JSON/JSONC formatting.
Most editors have an _autoformat on save_ feature — consider enabling
it for [Ruff][ruff-editors] and [Biome][biome-editors].

[pre-commit]: https://pre-commit.com/
[pre-commit.ci]: https://pre-commit.ci/
[Ruff]: https://docs.astral.sh/ruff/
[Biome]: https://biomejs.dev/
[ruff-editors]: https://docs.astral.sh/ruff/integrations/
[biome-editors]: https://biomejs.dev/guides/integrate-in-editor/

(writing-tests)=

## Writing tests

CellRank uses [pytest] for automated testing.
Please write {doc}`scanpy:dev/testing` for every function added to the package.

Most IDEs integrate with pytest and provide a GUI to run tests.
Point yours to one of the environments returned by

```bash
hatch env create hatch-test  # create test environments
hatch env find hatch-test    # list environment paths
```

Alternatively, you can run all tests from the command line:

:::::{tab-set}
::::{tab-item} Hatch
:sync: hatch

```bash
hatch test  # test with the highest supported Python version
# or
hatch test --all  # test with all supported Python versions
```

::::

::::{tab-item} uv
:sync: uv

```bash
uv run pytest
```

::::

::::{tab-item} Pip
:sync: pip

```bash
source .venv/bin/activate
pytest
```

::::
:::::

[pytest]: https://docs.pytest.org/

### Continuous integration

Continuous integration via GitHub Actions will automatically run the tests on all
pull requests, testing against the minimum and maximum supported Python version.

Additionally, there's a CI job that tests against pre-releases of all
dependencies (if there are any). This detects incompatibilities early and
gives you time to fix the issue or reach out to the dependency's developers
before the package is released to a wider audience.

The CI job is defined in `.github/workflows/test.yaml`,
however the single point of truth for CI jobs is the Hatch test matrix
defined in `pyproject.toml`. Local testing via hatch and remote CI testing
use the same Python versions and environments.

## Writing documentation

Please write documentation for new or changed features and use-cases.
This project uses [Sphinx] with the following features:

- [MyST] allows writing documentation in Markdown / Markedly Structured Text
- [NumPy-style docstrings][numpydoc] (through the [Napoleon][napoleon] extension)
- Jupyter notebooks as tutorials through [myst-nb]
  (see [Tutorials with myst-nb](#tutorials-with-myst-nb-and-jupyter-notebooks))
- [sphinx-autodoc-typehints] to automatically reference annotated input and output types
- Citations (like {cite:p}`lange:22`) via [sphinxcontrib-bibtex]

See scanpy's {doc}`scanpy:dev/documentation` for more information on how to write your own.

[Sphinx]: https://www.sphinx-doc.org/en/master/
[MyST]: https://myst-parser.readthedocs.io/en/latest/intro.html
[myst-nb]: https://myst-nb.readthedocs.io/en/latest/
[napoleon]: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
[numpydoc]: https://numpydoc.readthedocs.io/en/latest/format.html
[sphinx-autodoc-typehints]: https://github.com/tox-dev/sphinx-autodoc-typehints
[sphinxcontrib-bibtex]: https://sphinxcontrib-bibtex.readthedocs.io/

### Tutorials with myst-nb and Jupyter notebooks

The documentation renders Jupyter notebooks stored in the `docs/notebooks`
directory using [myst-nb]. Currently, only notebooks in `.ipynb` format are
supported that will be included with both their input and output cells.
It is your responsibility to update and re-run the notebook whenever necessary.

#### Hints

- If you refer to objects from other packages, add an entry to `intersphinx_mapping` in `docs/conf.py`.
  Only then can Sphinx automatically create a link to the external documentation.
- If building the documentation fails because of a missing link that is outside your control,
  you can add an entry to the `nitpick_ignore` list in `docs/conf.py`.

(building-the-docs-locally)=

### Building the docs locally

:::::{tab-set}
::::{tab-item} Hatch
:sync: hatch

```bash
hatch run docs:build
hatch run docs:open
```

For live preview with auto-rebuild on save:

```bash
hatch run docs:auto
```

::::

::::{tab-item} uv
:sync: uv

```bash
uv sync --group docs
source .venv/bin/activate
sphinx-build -M html docs docs/_build
open docs/_build/html/index.html  # macOS; use xdg-open on Linux
```

For live preview with auto-rebuild on save:

```bash
sphinx-autobuild docs docs/_build/html --open-browser --watch README.md
```

::::

::::{tab-item} Pip
:sync: pip

```bash
source .venv/bin/activate
cd docs
sphinx-build -M html . _build
open _build/html/index.html  # macOS; use xdg-open on Linux
```

::::
:::::

## Publishing a release

CellRank uses [hatch-vcs] for versioning — the version number is derived automatically from git tags.
To make a new release:

1. Navigate to the [Releases](https://github.com/theislab/cellrank/releases) page on GitHub.
2. Click "Draft a new release".
3. Specify a tag name of the form `vX.Y.Z`, adhering to [Semantic Versioning][semver].
4. Create the release — this triggers a GitHub workflow that builds and publishes to [PyPI].

After the PyPI release, the [conda-forge](https://conda-forge.org/) bot will automatically open a
PR to update the [CellRank feedstock](https://github.com/conda-forge/cellrank-feedstock).
Dependencies are synced from PyPI via Grayskull and the PR is auto-merged once CI passes.

[hatch-vcs]: https://pypi.org/project/hatch-vcs/
[semver]: https://semver.org/
[PyPI]: https://pypi.org/

## PR workflow

1. Fork the repository and create a feature branch from `main`.
2. Implement your changes, including tests and documentation.
3. Ensure pre-commit passes: `pre-commit run --all-files`.
4. Open a pull request against `main`.
5. CI will run tests automatically — check the results.
6. Address review feedback, then we'll merge!

Feel free to open a [draft PR](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests#draft-pull-requests)
early to get feedback on your approach.
If you have questions, open an [issue](https://github.com/theislab/cellrank/issues/new/choose)
or reach out via [Discourse](https://discourse.scverse.org/c/ecosystem/cellrank/).
