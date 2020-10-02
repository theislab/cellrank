# -*- coding: utf-8 -*-
"""Attribute, property and function names for estimators."""
from cellrank.tl._constants import PrettyEnum

META_KEY = "__prop_metadata__"


# attributes
class A(PrettyEnum):
    """Attribute names of an estimator."""

    EIG = "_eig"
    SCHUR = "_schur"
    SCHUR_MAT = "_schur_mat"
    MACRO = "_macrostates"
    MACRO_MEMBER = "_macrostates_memberships"
    MACRO_COLORS = "_macrostates_colors"
    TERM = "_terminal_states"
    TERM_PROBS = "_terminal_states_probabilities"
    TERM_COLORS = "_terminal_colors"
    TERM_ABS_PROBS = "_terminal_abs_probabilities"  # Lineage object
    ABS_PROBS = "_absorption_probabilities"
    DIFF_POT = "_diff_potential"
    COARSE_T = "_coarse_T"
    COARSE_INIT_D = "_coarse_init_dist"
    COARSE_STAT_D = "_coarse_stat_dist"
    LIN_ABS_TIMES = "_lin_abs_times"
    LIN_DRIVERS = "_lin_drivers"


# properties
class P(PrettyEnum):
    """Property names of an estimator."""

    NO_PROPERTY = "NO_PROPERTY"
    EMPTY = "EMPTY"
    EIG = "eigendecomposition"
    SCHUR = "schur"
    SCHUR_MAT = "schur_matrix"
    MACRO = "macrostates"
    MACRO_MEMBER = "macrostates_memberships"
    TERM = "terminal_states"
    TERM_PROBS = "terminal_states_probabilities"
    ABS_PROBS = "absorption_probabilities"
    DIFF_POT = "diff_potential"
    COARSE_T = "coarse_T"
    COARSE_INIT_D = "coarse_initial_distribution"
    COARSE_STAT_D = "coarse_stationary_distribution"
    LIN_ABS_TIMES = "lineage_absorption_times"
    LIN_DRIVERS = "lineage_drivers"


# functions
class F(PrettyEnum):
    """Function formats."""

    NO_FUNC = "NO_FUNC"
    GENERAL = "{}"
    COMPUTE = "compute_{}"
    PLOT = "plot_{}"
