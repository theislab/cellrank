About CellRank
==============
CellRank is a unified solution for the probabilistic description of cellular dynamics, encompassing various input data
modalities and analysis scenarios through one consistent application user interface (:doc:`API <../api/index>`).
If you find CellRank useful for your research, please check out :doc:`citing CellRank <cite>`.

Design principles
-----------------
Our framework is based on three key principles: `Robustness`_, `Modularity`_ and `Scalability`_.

Robustness
~~~~~~~~~~
Fate restriction is a gradual, noisy process requiring probabilistic treatment. Therefore, we use Markov chains to
describe stochastic fate transitions, where each cell represents one state in the Markov chain. Markov chains describe
memoryless transitions between system states through a probabilistic framework. Transition probabilities are summarized
in a transition matrix :math:`T`, where :math:`T_{ij}`, describes the probability of transitioning from state :math:`i`
to state :math:`j` in one step. In our context, each state corresponds to a cell. Markov chains are established tools in
single-cell genomics and form the basis of many successful pseudotime approaches :cite:`haghverdi:16,setty:19`.

By using Markov chains, we assume that cellular state transitions occur gradually and without memory. The former
assumption implies that cells change their molecular state in small steps with many intermediate states which are
captured in the data. This is a reasonable assumption for most biological systems. The latter assumption implies that a
state transition depends only on the current molecular state and not on the history of states. This assumption is valid
as CellRank describes average cellular dynamics, rather than any individual cell. Both assumptions form the basis of
many of the previous successful trajectory inference approaches :cite:`haghverdi:16,setty:19,wolf:19`.

Modularity
~~~~~~~~~~
A typical CellRank workflow consists of two steps: **(i)** estimating cell-cell transition probabilities to set up a
Markov transition matrix :math:`T`, and **(ii)** analyzing it using various tools to derive biological insights.
Decoupling these two steps yields a powerful and flexible modeling framework as many analysis steps are independent
of the construction of the transition matrix. For example, whether we use RNA velocity or a pseudotime to derive
directed transition probabilities does not change how initial and terminal states are inferred or fate probabilities
estimated. The general structure of the framework, corresponding to steps **(i)** and **(ii)**, is given by:

* :mod:`~cellrank.kernels` that take multi-view single cell input data  and estimate a matrix of cell-cell
  transition probabilities :math:`T`. Row :math:`i` in matrix :math:`T` contains the transition probabilities from cell
  :math:`i` towards putative descendants. Therefore, all entries in the matrix are between 0 and 1, and rows sum to one.
* :mod:`~cellrank.estimators` that take a cell-cell transition matrix :math:`T` computed using any kernel and
  apply concepts from the theory of Markov chains to identify initial, terminal, and intermediate macrostates
  and compute fate probabilities.

Our main (and recommended!) estimator is based on *Generalized Perron Cluster Cluster Analysis* (GPCCA)
:cite:`reuter:18,reuter:19`, a method originally developed to study molecular dynamics. CellRank uses a
robust implementation of GPCCA through the `pyGPCCA`_ package. Please don't forget to cite both CellRank and GPCCA when
using the :class:`~cellrank.estimators.GPCCA` estimator, see :doc:`citing CellRank <cite>`.

We use fate probabilities to visualize trajectory-specific gene expression trends, infer putative driver genes,
arrange cells in a :func:`circular embedding <cellrank.pl.circular_projection>` :cite:`velten:17`, visualize
:func:`cascades of gene activation <cellrank.pl.heatmap>` along a trajectory, and
:func:`cluster expression trends <cellrank.pl.cluster_trends>`. See our :doc:`tutorials <../notebooks/tutorials/index>` to learn more.

Scalability
~~~~~~~~~~~
All CellRank kernels yield sparse transition matrices :math:`T`. Further, the :class:`cellrank.estimators.GPCCA`
estimator exploits sparsity in all major computations. Sparsity allows CellRank to scale to millions of cells.

For example, when computing :meth:`fate probabilities <cellrank.estimators.GPCCA.compute_fate_probabilities>`, we transform the matrix
inversion problem into a set of linear systems, which we solve in parallel using the sparsity-optimized `GMRES`_ algorithm, implemented
efficiently in `PETSc`_. We use similar tricks to infer macrostates of cellular dyanmics via sparsity-optimized partial real Schur
decompositions (implemented under the hood via `pyGPCCA`_ and `SLEPc`_).

Why is it called "CellRank"?
----------------------------
CellRank **does not** rank cells, we gave the package this name because just like Google's original `PageRank`_
algorithm, it works with Markov chains to aggregate relationships between individual objects (cells vs. websites)
to learn about more global properties of the underlying dynamics (initial & terminal states and fate probabilities vs.
website relevance).

.. _PageRank: https://en.wikipedia.org/wiki/PageRank
.. _pyGPCCA: https://github.com/msmdev/pyGPCCA
.. _GMRES: https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
.. _PETSc: https://petsc.org/release/
.. _SLEPc: https://slepc.upv.es/
