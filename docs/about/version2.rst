Moving from CellRank v1 to v2
=============================
The original CellRank version (v1) (CITE) was a framework to analyze cellular dynamics based on RNA velocity (CITE) and
gene expression similarity. Based on a Markov chain (REF) formulation, it combined these two data modalities in a
high-dimensional space and used the GPCCA algorithm (REF, CITE) to compute initial and terminal states, fate probabilities, and driver
genes.

With version 2, we are generalizing CellRank beyond RNA velocity data and turn it into a general framework for single-cell
fate mapping based on various data modalities (see about CellRank (REF)). This required us to make substantial changes to CellRank's API, many of
which are **backward-compatability breaking**. CellRank 2 can do everything that version 1 could do (and much more);
however, it does it differently. We outline the most important changes here and refer to the release notes (REF) for a detailed
account.

Important changes in version 2
------------------------------
* Deprecation of high level `.tl` functions, including `cellrank.tl.terminal_states` and `cellrank.tl.lineages`: we realized
  that this mode of interacting with CellRank is not flexible enough to accomodate various data modalities and analysis paradigms.
  Thus, we switched to a more modular structure: :mod:`~cellrank.kernels` to compute cell-cell transition matrices,
  and :mod:`~cellrank.estimators` to compute initial and terminal states, fate probabilities, and more. See the tutorial (REF).
* Introduction of a :meth:`~cellrank.estimators.GPCCA.fit` :meth:`~cellrank.estimators.GPCCA.predict` workflow for estimators:
  given that we removed the old `.tl` high-level functions, we wanted to make it easier to interact with estimators. Thus,
  very estimator now has a `.fit` method, which computes macrostates, and a `.predict` method, which classifies some of these as
  terminal states. The new `.fit` `.predict` workflow complements our fully-flexible low-level mode of interacting with estimators.
  See the tutorial: (REF) (CITE)
* Renaming of `absorption_probabilities` to `fate_probabilities` everwhere, for example in :meth:`cellrank.estimators.GPCCA.compute_fate_probabilities`:
  while we still compute [absorption probabilities]() on the Markov chain under the hood, we realized that the term is somewhat technical and decided
  to replace it with the more intuitive `fate probabilities`.
* Removal of the `external API`: anyone wishing to contribute to CellRank can do this now directly via :mod:`~cellrank.kernels` and
  :mod:`~cellrank.estimators`. We welcome any contribution to CellRank, see our contribution guide (REF), and feel free to
  get in touch via an issue or email (REF).
* Replacement of the old `WOTKernel` with a new :class:`~ellrank.kernels.TransportMapKernel`: this is CellRank's interface
  with moscot (REF), enabling us to analyze large-scale time-course studies with additional spatial or lineage readout (CITE). In addition,
  the :class:`~ellrank.kernels.TransportMapKernel` interfaces with Waddington-OT (REF).

There are many more changes and improvements in CellRank 2. For example, the computation of fate probabilities is 30x faster compared
to version 1, we fixed many bugs, and improved and extended our documentation and tutorials (REF).
