# Installation

{mod}`cellrank` requires Python ≥ 3.12.
If you encounter any problems, feel free to open an [issue].

## pip / uv (recommended)

Install CellRank from [PyPI](https://pypi.org/project/cellrank/):

```bash
pip install cellrank
```

Or, using [uv](https://docs.astral.sh/uv/), a fast, modern pip replacement:

```bash
uv pip install cellrank
```

### Optional extras

CellRank provides optional dependency groups for common workflows:

```bash
pip install "cellrank[gam]"     # pyGAM for gene trends (cellrank.models.GAM)
pip install "cellrank[moscot]"  # moscot for optimal transport (RealTimeKernel)
pip install "cellrank[plot]"    # additional plotting dependencies
```

### PETSc and SLEPc (faster large-scale solvers)

For large datasets, CellRank uses [PETSc] and [SLEPc] (via [pyGPCCA]) to compute
macrostates and fate probabilities faster and with less memory than the default
[SciPy]-based eigensolver. Install them together with CellRank via the `petsc`
extra:

```bash
pip install "cellrank[petsc]"
```

`petsc4py` and `slepc4py` ship no pre-built wheels, so pip compiles PETSc and
SLEPc from source. CellRank only needs a serial build, which requires a **C
compiler** — no Fortran or MPI. On systems without a Fortran compiler or system
MPI (e.g. a stock macOS or minimal Linux), pass configuration options so the
build stays C-only:

```bash
PETSC_CONFIGURE_OPTIONS="--with-fc=0 --with-mpi=0 --with-debugging=0 --with-shared-libraries=1" \
    pip install "cellrank[petsc]"
```

Do **not** set `SLEPC_CONFIGURE_OPTIONS` — SLEPc inherits its configuration from
the PETSc build. The build is a one-time compilation (a few minutes) and is then
cached. If you would rather not compile from source, install PETSc and SLEPc via
conda instead (see below). If you specifically need a distributed **MPI** build,
install a system MPI and add `mpi4py` alongside the extra.

## conda / mamba / pixi

CellRank is available on [conda-forge](https://anaconda.org/conda-forge/cellrank).
This installation method is convenient if you want **pre-built** [PETSc] and
[SLEPc] — libraries for large-scale linear algebra that CellRank uses when
computing macrostates or fate probabilities on large datasets — without
compiling them from source (see the pip instructions above for the
compile-from-source alternative):

```bash
mamba install -c conda-forge cellrank
```

Or, using [pixi](https://pixi.sh/), a fast, modern conda replacement:

```bash
pixi add cellrank
```

To also install PETSc and SLEPc as pre-built conda packages (no compilation):

```bash
mamba install -c conda-forge cellrank petsc4py slepc4py
```

```{note}
**Windows users:** [PETSc] and [SLEPc] are not built for Windows on conda-forge
(only Linux and macOS). Because the conda-forge `cellrank` package depends on
`pygpcca`, which in turn requires `petsc`, `mamba install -c conda-forge cellrank`
cannot be solved on Windows and fails with an error such as
`petsc [...] does not exist (perhaps a missing channel)`.

On Windows, install CellRank with `pip install cellrank` instead. PETSc and SLEPc
are optional — without them CellRank falls back to a slower [SciPy]-based
eigensolver for computing macrostates and fate probabilities, which is fine for
small to medium datasets. To use PETSc and SLEPc on Windows, run CellRank inside
[WSL] and install via conda there.
```

## Development version

To install the latest development version from GitHub:

```bash
pip install git+https://github.com/theislab/cellrank.git@main
```

See {doc}`contributing` for setting up a full development environment.

[issue]: https://github.com/theislab/cellrank/issues/new
[PETSc]: https://petsc.org/
[SLEPc]: https://slepc.upv.es/
[pyGPCCA]: https://pygpcca.readthedocs.io/
[SciPy]: https://scipy.org/
[WSL]: https://learn.microsoft.com/windows/wsl/install
