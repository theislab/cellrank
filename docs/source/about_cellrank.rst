About CellRank
===============
CellRank is a unified solution for the probabilistic description of cellular dynamics, encompassing various input data
modalities and analysis scenarios through one consistent application user interface (API). Our framework is based on
three key principles: robustness, modularity and sparsity.

Robustness
----------
Fate restriction is a gradual, noisy process requiring probabilistic treatment. Therefore, we use Markov chains to
describe stochastic fate transitions, where each cell represents one state in the Markov chain. Markov chains describe
memoryless transitions between system states through a probabilistic framework. Transition probabilities are summarized
in a  transition matrix $T$, where $T_{jk}$, describes the probability of transitioning from state $j$ to state $k$ in
one step. In our context, each state corresponds to a cell. Markov chains are established tools in single-cell genomics
and form the basis of many successful pseudotime approaches :cite:`haghverdi:16,setty:19`.

By using Markov chains, we assume that cellular state transitions occur gradually and without memory. The former
assumption implies that cells change their molecular state in small steps with many intermediate states which are
captured in the data. This is a reasonable assumption for most biological systems. The latter assumption implies that a
state transition depends only on the current molecular state and not on the history of states. This assumption is valid
as CellRank describes average cellular dynamics, rather than any individual cell. Both assumptions form the basis of
many of the previous successful TI approaches :cite:`haghverdi:16,setty:19,wolf:19`.


Modularity
-----------
A typical CellRank workflow consists of two steps: (i) estimating cell-cell transition probabilities to set up a Markov
transition matrix $T$, and (ii) analyzing it using various tools to derive biological insights. Decoupling these two
steps yields a powerful and flexible modeling framework as many analysis steps are independent of the construction of
the transition matrix. For example, whether we use RNA velocity or a pseudotime to derive directed transition
probabilities does not change how initial and terminal states are inferred or fate probabilities estimated. The general
structure of the framework, corresponding to steps (i) and (ii), is given by:

* :ref:`Kernels <kernels>` that take multi-view single cell input data  and estimate a matrix of cell-cell transition
  probabilities $T$. Row $i$ in matrix $T$ contains the transition probabilities from cell $i$ towards putative
  descendants. Therefore, all entries in the matrix are between 0 and 1, and rows sum to one.
* :ref:`Estimators <estimators>` that take a cell-cell transition matrix $T$ computed using any kernel and apply
  concepts from the theory of Markov chains to identify initial, terminal, and intermediate `macrostates`_ and compute
  `fate probabilities`_.

Our main (and recommended!) estimator is based on *Generalized Perron Cluster Cluster Analysis* (GPCCA)
:cite:`reuter:18,reuter:19`, a method originally developed to study conformational protein dynamics. CellRank uses a
robust implementation of GPCCA through the `pyGPCCA`_ package.
Please don't forget to cite both CellRank and GPCCA when using the :class:`cellrank.estimators.GPCCA` estimator,
see :doc:`How to cite CellRank <citing_cellrank>`.

We use fate probabilities to visualize trajectory-specific `gene expression trends`_, infer putative `driver genes`_,
arrange cells in a `circular embedding`_ :cite:`velten:17`,  visualize `cascades of gene activation`_ along a
trajectory, and `cluster expression trends`_.

Sparsity
--------
All CellRank kernels yield sparse transition matrices $T$. Further, the :class:`cellrank.estimators.GPCCA` estimator
exploits sparsity in all major computations. Sparsity allows CellRank to scale to large datasets.


Why is it called "CellRank"?
----------------------------
CellRank **does not** rank cells, we gave the package this name because just like Google's original `PageRank`_
algorithm, it works with Markov chains to aggregate relationships between individual objects (cells vs. websites)
to learn about more global properties of the underlying dynamics (initial & terminal states and fate probabilities vs.
website relevance).


.. _PageRank: https://en.wikipedia.org/wiki/PageRank#cite_note-1
.. _pyGPCCA: https://pygpcca.readthedocs.io/

.. _macrostates: :doc:`notebooks/tutorials/initial_terminal_states`
.. _fate probabilities: :doc:`notebooks/tutorials/fate_probabilities`
.. _driver genes: :doc:`notebooks/tutorials/fate_probabilities`
.. _gene expression trends: :doc:`notebooks/tutorials/gene_trends`
.. _circular embedding: :func:`cellrank.pl.circular_projection`
.. _cascades of gene activation: :func:`cellrank.pl.heatmap`
.. _cluster expression trends: :func:`cellrank.pl.cluster_trends`
