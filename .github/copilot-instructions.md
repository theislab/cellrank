# Copilot Instructions for CellRank

## What CellRank Does

CellRank analyzes **cellular dynamics** from single-cell data by modeling cells
as states in a Markov chain. It estimates transition probabilities between cells,
identifies initial and terminal cell states, and computes fate probabilities —
the likelihood that each cell will reach a given terminal state.

**Two-layer architecture**:
1. **Kernels** compute cell-cell transition matrices from biological signals
   (RNA velocity, pseudotime, real time + optimal transport, etc.)
2. **Estimators** analyze these transition matrices to find macrostates, terminal
   states, and fate probabilities using spectral methods (GPCCA)

## Important Notes

- Activate environment before running code: `source .venv/bin/activate`
- **Run tests with tox** (not pytest directly): `tox -e py312-noslepc`
  or specific test files: `tox -e py312-noslepc -- tests/test_kernels.py -v`
  (will migrate to `hatch test` — see `.github/prompts/PLAN_modernize.md`)

## Module Guide

### `kernels/` — Transition Matrix Construction

Each kernel turns a biological signal into a row-stochastic transition matrix $T$
where $T_{ij}$ is the probability of cell $i$ transitioning to cell $j$.

Kernels are either **bidirectional** or **unidirectional**:
- **Bidirectional** kernels support both forward and backward analysis.
  `~kernel` flips the direction, enabling study of both fate commitment
  (forward) and cellular origin (backward).
- **Unidirectional** kernels only compute transitions in one direction
  (e.g., `RealTimeKernel` only moves forward in real time).

| Kernel | Signal | Direction | Key deps |
|--------|--------|-----------|----------|
| `VelocityKernel` | Cell-level vector fields (RNA velocity, flow matching, Wasserstein gradient flows, etc.) | Bidirectional | scvelo; jax (stochastic mode) |
| `PseudotimeKernel` | Pseudotime ordering | Bidirectional | — |
| `CytoTRACEKernel` | Gene count → plasticity score | Bidirectional | — |
| `RealTimeKernel` | OT couplings across timepoints | Unidirectional | moscot or wot |
| `ConnectivityKernel` | kNN graph similarity | Unidirectional (symmetric) | — |
| `PrecomputedKernel` | Any external transition matrix | Unidirectional | — |

**Kernel composition** — kernels combine via arithmetic operators:
```python
vk = VelocityKernel(adata).compute_transition_matrix()
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined = 0.8 * vk + 0.2 * ck  # KernelAdd(KernelMul(...), KernelMul(...))
```
This builds an expression tree. `+` normalizes weights to sum to 1, `*` does
element-wise product.

**Key base classes** (in `_base_kernel.py`):
- `Kernel` — abstract root with `transition_matrix` property, `write_to_adata()`,
  plotting methods
- `UnidirectionalKernel` / `BidirectionalKernel` — leaf kernel bases
- `KernelExpression` / `KernelAdd` / `KernelMul` — composite expression nodes
- `Constant` — wraps a scalar weight

**`RealTimeKernel`** is the most complex kernel — it assembles couplings from
moscot (`from_moscot()`) or WOT (`from_wot()`) into a global transition matrix
via `_restitch_couplings()`, with configurable self-transition blocks.

### `estimators/` — Markov Chain Analysis

| Estimator | Method | Use case |
|-----------|--------|----------|
| `GPCCA` | Generalized Perron Cluster Cluster Analysis | **Recommended** — spectral soft clustering |
| `CFLARE` | Clustering/Filtering of Left/Right Eigenvectors | Legacy method |

**GPCCA workflow** (the standard path):
```python
g = cr.estimators.GPCCA(kernel)
g.compute_schur(n_components=10)
g.compute_macrostates(n_states=5)
g.predict_terminal_states()
# OR: g.set_terminal_states(states=["Alpha", "Beta"])
g.compute_fate_probabilities()
g.compute_lineage_drivers()
```

Many users set terminal states manually via `set_terminal_states()` based on
prior biological knowledge, rather than relying on automatic prediction.

**Mixin architecture** — `GPCCA` inherits from 4 mixins in `estimators/mixins/`:
- `SchurMixin` (`decomposition/_schur.py`) — Schur decomposition
- `EigenMixin` (`decomposition/_eigen.py`) — standard eigendecomposition
- `FateProbsMixin` (`_fate_probabilities.py`) — absorption probabilities
- `LineageDriversMixin` (`_lineage_drivers.py`) — gene-fate correlations

### `models/` — Gene Trend Fitting

Models fit smooth curves to gene expression along lineages for visualization.
Pipeline: `prepare()` → `fit()` → `predict()` → `confidence_interval()`.

| Model | Backend | Optional dep |
|-------|---------|--------------|
| `GAM` | pygam (Python) | — (hard dep) |
| `GAMR` | R's mgcv via rpy2 | rpy2 + R |
| `SKLearnModel` | Any sklearn estimator | — |

### `pl/` — Plotting

Six main plotting functions: `aggregate_fate_probabilities`, `circular_projection`,
`cluster_trends`, `gene_trend`, `heatmap`, `log_odds`.

### `_utils/` — Core Utilities

- **`Lineage`** (`_lineage.py`) — NumPy ndarray subclass for fate probabilities
  (cells × terminal states). Has named columns, colors, distance methods.
- **`_linear_solver.py`** — iterative solvers (GMRES, BiCGSTAB) + optional PETSc/SLEPc
- **`_parallelize.py`** — joblib-based parallelization helper
- **`_key.py`** — standardized AnnData key naming (e.g., `terminal_states_fwd`)

### `logging/` — Custom Logger

Legacy scanpy-style logger (`_RootLogger`) with custom `HINT` level, timing
support, and `settings.verbosity` integration. Used as:
```python
from cellrank import logging as logg

logg.info("Computing transition matrix", time=start)
```
35+ files use this pattern — don't refactor it.

### `datasets.py` — Example Datasets

Downloads from Figshare: `pancreas` (endocrinogenesis), `lung` (regeneration),
`reprogramming_morris`, `reprogramming_schiebinger`, `zebrafish`, `bone_marrow`.

## pyGPCCA — Sibling Repo

CellRank's GPCCA estimator delegates to the
[pyGPCCA](https://github.com/msmdev/pyGPCCA) library (installed as `pygpcca`,
local checkout at `../pyGPCCA/`). It implements the GPCCA algorithm: sorted real
Schur decomposition followed by rotation optimization to identify metastable
(macro-) states as fuzzy clusters.

## Testing

**Test structure**: flat layout in `tests/`, one file per major component.

| File | Tests |
|------|-------|
| `test_kernels.py` | All kernel types, composition, I/O |
| `test_gpcca.py` | GPCCA estimator workflow |
| `test_cflare.py` | CFLARE estimator |
| `test_lineage.py` | Lineage array operations |
| `test_model.py` | GAM/GAMR/sklearn model fitting |
| `test_plotting.py` | Figure comparison tests |
| `test_pipeline.py` | End-to-end workflows |
| `test_linear_solver.py` | Sparse solvers |
| `test_lineage_drivers.py` | Gene-fate correlations |

**Test data**: pre-computed h5ad files in `tests/_ground_truth_adatas/`
(50, 100, 200 cells). Fixtures in `conftest.py` pre-build CFLARE/GPCCA
estimators with `VelocityKernel + ConnectivityKernel`.

**Running tests**:
```bash
source .venv/bin/activate
tox -e py312-noslepc                          # standard tests
tox -e py312-slepc                            # PETSc/SLEPc tests (needs conda)
python -m pytest tests/test_kernels.py -v -x  # quick iteration
```

**109 skips expected**: PETSc/SLEPc and R tests skip when deps not installed.

## Key Patterns

- **AnnData-centric**: kernels read from / write to `adata.obsp`, `adata.obs`,
  `adata.uns`. Estimators store results in the same AnnData.
- **`write_to_adata()`**: persists transition matrix + params to AnnData
- **`from_adata()` classmethod**: reconstructs kernel/estimator from saved AnnData
- **Shadow AnnData**: estimators maintain an internal AnnData copy for
  serialization; `to_adata()` exports it
- **Bidirectional kernels**: `~kernel` flips direction (forward ↔ backward)
- **Logging**: Use `from cellrank import logging as logg` (not stdlib logging)
- **Key naming**: Use helpers in `_utils/_key.py` for AnnData keys

## Optional Dependencies

| Package | Used by | Purpose |
|---------|---------|---------|
| `scvelo` | VelocityKernel, CytoTRACEKernel | RNA velocity, moments |
| `moscot` | RealTimeKernel.from_moscot() | OT couplings (JAX/OTT backend) |
| `jax` / `jaxlib` | VelocityKernel (stochastic) | Hessian computation |
| `petsc4py` / `slepc4py` | SchurMixin, FateProbsMixin | Krylov Schur, large linear systems |
| `rpy2` + R (`mgcv`) | GAMR model | R-based GAMs |
| `adjustText` | lineage drivers plot | Label placement |
| `wot` | RealTimeKernel.from_wot() | Waddington-OT format |

Currently no pip extras for these — install manually.

## Related Projects

> Use `github_repo` tool to search these when relevant topics come up.

- **[moscot](https://github.com/theislab/moscot)** — Optimal transport for
  single-cell (JAX/OTT). Computes couplings used by `RealTimeKernel`.
- **[pyGPCCA](https://github.com/msmdev/pyGPCCA)** — GPCCA algorithm
  implementation. Core dependency for macrostate identification.
- **[scvelo](https://github.com/theislab/scvelo)** — RNA velocity.
  Used by `VelocityKernel` for velocity computation.
- **[scanpy](https://github.com/scverse/scanpy)** — Single-cell analysis
  toolkit. CellRank's foundation (kNN graph, preprocessing).
- **[scverse cookiecutter](https://github.com/scverse/cookiecutter-scverse)** —
  Template this package is migrating toward.
