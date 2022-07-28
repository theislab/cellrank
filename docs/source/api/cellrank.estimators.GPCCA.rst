.. taken from https://github.com/scverse/scvi-tools/blob/master/docs/_templates/autosummary/class.rst

cellrank.estimators.GPCCA
=========================

.. currentmodule:: cellrank.estimators
.. autoclass:: GPCCA




Attributes table
~~~~~~~~~~~~~~~~
.. autosummary::


    ~cellrank.estimators.GPCCA.absorption_probabilities


    ~cellrank.estimators.GPCCA.absorption_times


    ~cellrank.estimators.GPCCA.adata


    ~cellrank.estimators.GPCCA.backward


    ~cellrank.estimators.GPCCA.coarse_T


    ~cellrank.estimators.GPCCA.coarse_initial_distribution


    ~cellrank.estimators.GPCCA.coarse_stationary_distribution


    ~cellrank.estimators.GPCCA.eigendecomposition


    ~cellrank.estimators.GPCCA.kernel


    ~cellrank.estimators.GPCCA.lineage_drivers


    ~cellrank.estimators.GPCCA.macrostates


    ~cellrank.estimators.GPCCA.macrostates_memberships


    ~cellrank.estimators.GPCCA.params


    ~cellrank.estimators.GPCCA.priming_degree


    ~cellrank.estimators.GPCCA.schur_matrix


    ~cellrank.estimators.GPCCA.schur_vectors


    ~cellrank.estimators.GPCCA.shape


    ~cellrank.estimators.GPCCA.terminal_states


    ~cellrank.estimators.GPCCA.terminal_states_memberships


    ~cellrank.estimators.GPCCA.terminal_states_probabilities


    ~cellrank.estimators.GPCCA.transition_matrix








Methods table
~~~~~~~~~~~~~
.. autosummary::






    ~cellrank.estimators.GPCCA.compute_absorption_probabilities




    ~cellrank.estimators.GPCCA.compute_eigendecomposition




    ~cellrank.estimators.GPCCA.compute_lineage_drivers




    ~cellrank.estimators.GPCCA.compute_lineage_priming




    ~cellrank.estimators.GPCCA.compute_macrostates




    ~cellrank.estimators.GPCCA.compute_schur




    ~cellrank.estimators.GPCCA.compute_terminal_states




    ~cellrank.estimators.GPCCA.copy




    ~cellrank.estimators.GPCCA.fit




    ~cellrank.estimators.GPCCA.from_adata




    ~cellrank.estimators.GPCCA.plot_absorption_probabilities




    ~cellrank.estimators.GPCCA.plot_coarse_T




    ~cellrank.estimators.GPCCA.plot_lineage_drivers




    ~cellrank.estimators.GPCCA.plot_lineage_drivers_correlation




    ~cellrank.estimators.GPCCA.plot_macrostate_composition




    ~cellrank.estimators.GPCCA.plot_macrostates




    ~cellrank.estimators.GPCCA.plot_schur_matrix




    ~cellrank.estimators.GPCCA.plot_spectrum




    ~cellrank.estimators.GPCCA.plot_terminal_states




    ~cellrank.estimators.GPCCA.predict




    ~cellrank.estimators.GPCCA.read




    ~cellrank.estimators.GPCCA.rename_terminal_states




    ~cellrank.estimators.GPCCA.set_terminal_states




    ~cellrank.estimators.GPCCA.set_terminal_states_from_macrostates




    ~cellrank.estimators.GPCCA.to_adata




    ~cellrank.estimators.GPCCA.write









Attributes
~~~~~~~~~~


absorption\_probabilities
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.absorption_probabilities



absorption\_times
^^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.absorption_times



adata
^^^^^
.. autoattribute:: GPCCA.adata



backward
^^^^^^^^
.. autoattribute:: GPCCA.backward



coarse\_T
^^^^^^^^^
.. autoattribute:: GPCCA.coarse_T



coarse\_initial\_distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.coarse_initial_distribution



coarse\_stationary\_distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.coarse_stationary_distribution



eigendecomposition
^^^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.eigendecomposition



kernel
^^^^^^
.. autoattribute:: GPCCA.kernel



lineage\_drivers
^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.lineage_drivers



macrostates
^^^^^^^^^^^
.. autoattribute:: GPCCA.macrostates



macrostates\_memberships
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.macrostates_memberships



params
^^^^^^
.. autoattribute:: GPCCA.params



priming\_degree
^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.priming_degree



schur\_matrix
^^^^^^^^^^^^^
.. autoattribute:: GPCCA.schur_matrix



schur\_vectors
^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.schur_vectors



shape
^^^^^
.. autoattribute:: GPCCA.shape



terminal\_states
^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.terminal_states



terminal\_states\_memberships
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.terminal_states_memberships



terminal\_states\_probabilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.terminal_states_probabilities



transition\_matrix
^^^^^^^^^^^^^^^^^^
.. autoattribute:: GPCCA.transition_matrix









Methods
~~~~~~~






compute\_absorption\_probabilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.compute_absorption_probabilities





compute\_eigendecomposition
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.compute_eigendecomposition





compute\_lineage\_drivers
^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.compute_lineage_drivers





compute\_lineage\_priming
^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.compute_lineage_priming





compute\_macrostates
^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.compute_macrostates





compute\_schur
^^^^^^^^^^^^^^
.. automethod:: GPCCA.compute_schur





compute\_terminal\_states
^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.compute_terminal_states





copy
^^^^
.. automethod:: GPCCA.copy





fit
^^^
.. automethod:: GPCCA.fit





from\_adata
^^^^^^^^^^^
.. automethod:: GPCCA.from_adata





plot\_absorption\_probabilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.plot_absorption_probabilities





plot\_coarse\_T
^^^^^^^^^^^^^^^
.. automethod:: GPCCA.plot_coarse_T





plot\_lineage\_drivers
^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.plot_lineage_drivers





plot\_lineage\_drivers\_correlation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.plot_lineage_drivers_correlation





plot\_macrostate\_composition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.plot_macrostate_composition





plot\_macrostates
^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.plot_macrostates





plot\_schur\_matrix
^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.plot_schur_matrix





plot\_spectrum
^^^^^^^^^^^^^^
.. automethod:: GPCCA.plot_spectrum





plot\_terminal\_states
^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.plot_terminal_states





predict
^^^^^^^
.. automethod:: GPCCA.predict





read
^^^^
.. automethod:: GPCCA.read





rename\_terminal\_states
^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.rename_terminal_states





set\_terminal\_states
^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.set_terminal_states





set\_terminal\_states\_from\_macrostates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: GPCCA.set_terminal_states_from_macrostates





to\_adata
^^^^^^^^^
.. automethod:: GPCCA.to_adata





write
^^^^^
.. automethod:: GPCCA.write
