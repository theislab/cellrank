Moving to CellRank 2
====================
The original CellRank version (v1) :cite:`lange:22` was a framework to analyze cellular dynamics based on RNA velocity
:cite:`manno:18,bergen:20` and gene expression similarity. Based on a
`Markov chain <https://en.wikipedia.org/wiki/Markov_chain>`_ formulation, it combined these two data modalities
in a high-dimensional space and used the :class:`GPCCA algorithm <pygpcca.GPCCA>` :cite:`reuter:18,reuter:22`
to compute initial and terminal states, fate probabilities, and driver genes.

With version 2, we are generalizing CellRank beyond RNA velocity data and turn it into a general framework for
single-cell fate mapping based on various data modalities (see :doc:`index`). This required us to make substantial
changes to CellRank's API, many of which are **backward-compatibility breaking**. CellRank 2 can do everything that
version 1 could do (and much more); however, it does it differently.

.. important::
    We outline the most important changes here. For a more detailed account, please refer to the
    :doc:`../release_notes`.

Important changes in version 2
------------------------------
* Deprecation of high level ``cellrank.tl`` functions, including ``cellrank.tl.terminal_states`` and
  ``cellrank.tl.lineages``: we realized that this mode of interacting with CellRank is not flexible enough to
  accommodate various data modalities and analysis paradigms. Thus, we switched to a more modular structure:
  :mod:`cellrank.kernels` to compute cell-cell transition matrices, and :mod:`cellrank.estimators` to compute initial
  and terminal states, fate probabilities, and more. See :doc:`../notebooks/tutorials/general/100_getting_started`.
* Introduction of a :meth:`~cellrank.estimators.GPCCA.fit`/:meth:`~cellrank.estimators.GPCCA.predict` workflow for
  estimators: given that we removed the old ``cellrank.tl`` high-level functions, we wanted to make it easier to
  interact with estimators. Thus, every :mod:`estimator <cellrank.estimators>` now has a ``.fit`` method, which computes
  macrostates, and a ``.predict`` method, which classifies some of these as terminal states.
  The new :meth:`~cellrank.estimators.GPCCA.fit`/:meth:`~cellrank.estimators.GPCCA.predict` workflow complements our
  fully-flexible low-level mode of interacting with estimators.
  See :doc:`../notebooks/tutorials/estimators/600_initial_terminal`.
* Renaming of ``absorption_probabilities`` to :attr:`~cellrank.estimators.GPCCA.fate_probabilities` everywhere,
  for example in :meth:`~cellrank.estimators.GPCCA.compute_fate_probabilities`: while we still compute
  `absorption probabilities <https://en.wikipedia.org/wiki/Absorbing_Markov_chain>`_ on the Markov chain under the hood,
  we realized that the term is somewhat technical and decided to replace it with the more intuitive
  "fate probabilities". See :doc:`../notebooks/tutorials/estimators/700_fate_probabilities`.
* Removal of the ``cellrank.external``: anyone wishing to contribute to CellRank can do this now directly via
  :mod:`cellrank.kernels` and :mod:`cellrank.estimators`. We welcome any contribution to CellRank,
  see our :doc:`contribution guide <../contributing>`, and feel free to get in touch via an
  `issue <https://github.com/theislab/cellrank/issues/new/choose>`_ or `email <mailto:info@cellrank.org>`_.
* Replacement of the old ``cellrank.external.WOTKernel`` with a new :class:`cellrank.kernels.RealTimeKernel`: this is
  CellRank's interface with :mod:`moscot`, enabling us to analyze large-scale time-course studies with additional
  spatial or lineage readout :cite:`klein:23,lange:23`. In addition, the :class:`~cellrank.kernels.RealTimeKernel`
  interfaces with `Waddington-OT <https://broadinstitute.github.io/wot/>`_ :cite:`schiebinger:19`.

There are many more changes and improvements in CellRank 2. For example, the computation of fate probabilities is
**30x** faster compared to version 1, we fixed many bugs, and improved and extended our documentation and
:doc:`tutorials <../notebooks/tutorials/index>`.
