import textwrap
from typing import Any

from docrep import DocstringProcessor

__all__ = ["d", "inject_docs"]

_adata = """\
adata
    Annotated data object."""
_adata_ret = """\
:class:`anndata.AnnData`
    Annotated data object."""
_plotting = """\
figsize
    Size of the figure.
dpi
    Dots per inch.
save
    Filename where to save the plot."""
_n_jobs = """\
n_jobs
    Number of parallel jobs. If -1, use all available cores. If :obj:`None` or 1, the execution is sequential."""
_parallel = f"""\
show_progress_bar
    Whether to show a progress bar. Disabling it may slightly improve performance.
{_n_jobs}
backend
    Which backend to use for parallelization. See :class:`~joblib.Parallel` for valid options."""
_model = """\
model
    Model based on :class:`~cellrank.models.BaseModel` to fit.
    If a :class:`dict`, gene and lineage specific models can be specified. Use ``'*'`` to indicate
    all genes or lineages, for example ``{'gene_1': {'*': ...}, 'gene_2': {'lineage_1': ..., '*': ...}}``."""
_just_plots = """\
Nothing, just plots the figure. Optionally saves it based on ``save``."""
_plots_or_returns_models = """\
If ``return_models = False``, just plots the figure and optionally saves it based on ``save``.
Otherwise returns the fitted models as ``{'gene_1': {'lineage_1': <model_11>, ...}, ...}``.
Models which have failed will be instances of :class:`cellrank.models.FailedModel`."""
_backward = """\
backward
    Direction of the process."""
_eigen = """\
which
    How to sort the eigenvalues. Valid option are:

    - ``'LR'`` - the largest real part.
    - ``'LM'`` - the largest magnitude.
alpha
    Used to compute the *eigengap*. ``alpha`` is the weight given to the deviation of an eigenvalue from one."""
_n_cells = """\
n_cells
    Number of most likely cells from each macrostate to select."""
_fit = """\
n_lineages
    Number of lineages. If :obj:`None`, it will be determined automatically.
cluster_key
    Match computed states against pre-computed clusters to annotate the states.
    For this, provide a key from :attr:`anndata.AnnData.obs` where cluster labels have been computed.
keys
    Determines which states to use by passing their names. Furthermore, these states can be combined, e.g.,
    if the states are ``['Neuronal_1', 'Neuronal_1', 'OPC']``, then passing
    ``keys = ['Neuronal_1, Neuronal_2', 'OPC']`` means that the two neuronal states are treated as one."""
_density_correction = (
    "Optionally, we apply a density correction as described in :cite:`coifman:05`, "
    "where we use the implementation of :cite:`haghverdi:16`."
)
_time_range = """\
time_range
    Specify start and end times:

    - :class:`tuple` - it specifies the minimum and maximum pseudotime. Both values can be :obj:`None`,
      in which case the minimum is the earliest pseudotime and the maximum is automatically determined.
    - :class:`float` - it specifies the maximum pseudotime."""

_velocity_mode = """\
mode
    How to compute transition probabilities. Valid options are:

    - ``{m.DETERMINISTIC!r}`` - deterministic computation that doesn't propagate uncertainty.
    - ``{m.MONTE_CARLO!r}`` - Monte Carlo average of randomly sampled velocity vectors.
    - ``{m.STOCHASTIC!r}`` - second order approximation, only available when :mod:`jax` is installed."""
_velocity_backward_mode = """\
backward_mode
    Only matters if initialized as :attr:`backward = True <backward>`. Valid options are:

    - ``{b.TRANSPOSE!r}`` - compute transitions from neighboring cells :math:`j` to cell :math:`i`.
    - ``{b.NEGATE!r}`` - negate the velocity vector."""
_copy = """Return a copy of self."""
_initial = "initial"
_terminal = "terminal"
_model_callback = """\
callback
    Function which takes a :class:`~cellrank.models.BaseModel` and some keyword arguments
    for :meth:`~cellrank.models.BaseModel.prepare` and returns the prepared model.
    Can be specified in gene- and lineage-specific manner, similarly to the :attr:`model`."""
_genes = """\
genes
    Genes in :attr:`~anndata.AnnData.var_names`."""
_softmax_scale = """\
softmax_scale
    Scaling parameter for the softmax. If :obj:`None`, it will be estimated using ``1 / median(correlations)``.
    The idea behind this is to scale the softmax to counter everything tending to orthogonality in high dimensions."""
_time_mode = """\
mode
    Valid options are:

    - ``'embedding'`` - plot the embedding while coloring in continuous or categorical observations.
    - ``'time'`` - plot the pseudotime on x-axis and the probabilities/memberships on y-axis."""
_write_to_adata = """\
Updates the :attr:`adata` with the following fields:

- :attr:`obsp['{{key}}'] <anndata.AnnData.obsp>` - the transition matrix.
- :attr:`uns['{{key}}_params'] <anndata.AnnData.uns>` - parameters used for the calculation."""
_en_cutoff_p_thresh = """\
en_cutoff
    If ``cluster_key`` is given, this parameter determines when an approximate recurrent class will
    be labeled as *'Unknown'*, based on the entropy of the distribution of cells over transcriptomic clusters.
p_thresh
    If cell cycle scores were provided, a *Wilcoxon rank-sum test* is conducted to identify cell-cycle states.
    If the test returns a positive statistic and a p-value smaller than ``p_thresh``, a warning will be issued."""
_return_models = """\
return_models
    If :obj:`True`, return the fitted models for each gene in ``genes`` and lineage in ``lineages``."""
_basis = """\
basis
    Basis to use when ``mode = 'embedding'``. If :obj:`None`, use ``'umap'``."""
_velocity_scheme = """\
similarity
    Similarity measure between cells as described in :cite:`li:20`. Valid options are:

    - ``{s.CORRELATION!r}`` - :class:`~cellrank.kernels.utils.Correlation`.
    - ``{s.COSINE!r}`` - :class:`~cellrank.kernels.utils.Cosine`.
    - ``{s.DOT_PRODUCT!r}`` - :class:`~cellrank.kernels.utils.DotProduct`.

    Alternatively, any function can be passed as long as it follows the signature of
    :meth:`cellrank.kernels.utils.SimilarityABC.__call__`."""
_cond_num = """\
compute_cond_num
    Whether to compute condition number of the transition matrix. Note that this might be costly,
    since it does not use sparse implementation."""
_soft_scheme_fmt = """\
b
    The growth rate of generalized logistic function.{}
nu
    Affects near which asymptote maximum growth occurs.{}"""
_rw_ixs = """\
Can be specified as:

    - :class:`dict` - dictionary with 1 key in :attr:`~anndata.AnnData.obs` with values corresponding
      to either 1 or more clusters (if the column is categorical) or a :class:`tuple` specifying
      :math:`[min, max]` interval from which to select the indices.
    - :class:`~typing.Sequence` - sequence of cell ids in :attr:`~anndata.AnnData.obs_names`.
"""
_gene_symbols = """\
gene_symbols
    Key in :attr:`~anndata.AnnData.var` to use instead of :attr:`~anndata.AnnData.var_names`."""
_absorption_utils = """\
solver
    Solver to use for the linear problem. Options are ``'direct'``, ``'gmres'``, ``'lgmres'``, ``'bicgstab'`` or
    ``'gcrotmk'`` when ``use_petsc = False``.

    Information on the :mod:`scipy` iterative solvers can be found in :mod:`scipy.sparse.linalg` or for
    the ``petsc`` solvers `here <https://petsc.org/release/overview/linear_solve_table/>`__.
use_petsc
    Whether to use solvers from :mod:`petsc4py` or :mod:`scipy`. Recommended for large problems.
    If no installation is found, defaults to :func:`~scipy.sparse.linalg.gmres`.
n_jobs
    Number of parallel jobs to use when using an iterative solver.
backend
    Which backend to use for multiprocessing. See :class:`~joblib.Parallel` for valid options.
show_progress_bar
    Whether to show progress bar. Only used when ``solver != 'direct'``.
tol
    Convergence tolerance for the iterative solver. The default is fine for most cases, only consider
    decreasing this for severely ill-conditioned matrices.
preconditioner
    Preconditioner to use, only available when ``use_petsc = True``. For valid options, see
    `here <https://petsc.org/release/docs/manual/ksp/?highlight=pctype#preconditioners>`__.
    We recommend the ``'ilu'`` preconditioner for badly conditioned problems."""
which = """\
which
    Whether to compute initial or terminal states."""
allow_overlap = """\
allow_overlap
    Whether to allow overlapping names between initial and terminal states.
"""


def inject_docs(**kwargs: Any):
    def decorator(obj):
        obj.__doc__ = textwrap.dedent(obj.__doc__).format(**kwargs)
        return obj

    def decorator2(obj):
        obj.__doc__ = textwrap.dedent(kwargs["__doc__"])
        return obj

    if isinstance(kwargs.get("__doc__", None), str) and len(kwargs) == 1:
        return decorator2

    return decorator


d = DocstringProcessor(
    plotting=_plotting,
    n_jobs=_n_jobs,
    parallel=_parallel,
    model=_model,
    adata=_adata,
    adata_ret=_adata_ret,
    just_plots=_just_plots,
    backward=_backward,
    initial=_initial,
    terminal=_terminal,
    eigen=_eigen,
    initial_or_terminal=f"{_initial} or {_terminal}",
    n_cells=_n_cells,
    fit=_fit,
    copy=_copy,
    density_correction=_density_correction,
    time_range=_time_range,
    velocity_mode=_velocity_mode,
    velocity_backward_mode=_velocity_backward_mode,
    model_callback=_model_callback,
    genes=_genes,
    softmax_scale=_softmax_scale,
    time_mode=_time_mode,
    write_to_adata=_write_to_adata,
    en_cutoff_p_thresh=_en_cutoff_p_thresh,
    return_models=_return_models,
    plots_or_returns_models=_plots_or_returns_models,
    basis=_basis,
    velocity_scheme=_velocity_scheme,
    cond_num=_cond_num,
    soft_scheme=_soft_scheme_fmt.format("", "", ""),
    soft_scheme_kernel=_soft_scheme_fmt.format(*([" Only used when ``threshold_scheme = 'soft'``."] * 3)),
    rw_ixs=_rw_ixs,
    gene_symbols=_gene_symbols,
    absorption_utils=_absorption_utils,
    which=which,
    allow_overlap=allow_overlap,
)
