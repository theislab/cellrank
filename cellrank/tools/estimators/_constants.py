# -*- coding: utf-8 -*-
"""Attribute, property and function names for estimators."""
from cellrank.tools._constants import PrettyEnum

META_KEY = "__prop_metadata__"


# attributes
class A(PrettyEnum):
    """Attribute names of an estimator."""

    EIG = "_eig"
    SCHUR = "_schur"
    SCHUR_MAT = "_schur_mat"
    META = "_metastable_states"
    META_PROBS = "_metastable_states_probabilities"
    META_COLORS = "_metastable_states_colors"
    FIN = "_final_states"
    FIN_PROBS = "_final_states_probabilities"
    FIN_COLORS = "_final_colors"
    FIN_ABS_PROBS = "_final_abs_probabilities"  # the Lineage object
    ABS_RPOBS = "_absorption_probabilities"
    DIFF_POT = "_diff_potential"
    COARSE_T = "_coarse_T"
    COARSE_INIT_D = "_coarse_init_dist"
    COARSE_STAT_D = "_coarse_stat_dist"
    MEAN_ABS_TIME = "_mean_abs_time"
    VAR_ABS_TIME = "_var_abs_time"


# properties
class P(PrettyEnum):
    """Property names of an estimator."""

    NO_PROPERTY = "NO_PROPERTY"
    EMPTY = "EMPTY"
    EIG = "eigendecomposition"
    SCHUR = "schur"
    SCHUR_MAT = "schur_matrix"
    META = "metastable_states"
    META_PROBS = "metastable_states_probabilities"
    FIN = "final_states"
    FIN_PROBS = "final_states_probabilities"
    ABS_PROBS = "absorption_probabilities"
    DIFF_POT = "diff_potential"
    COARSE_T = "coarse_T"
    COARSE_INIT_D = "coarse_initial_distribution"
    COARSE_STAT_D = "coarse_stationary_distribution"
    MEAN_ABS_TIME = "absorption_time"
    VAR_ABS_TIME = "absorption_time_variance"


# functions
class F(PrettyEnum):
    """Function formats."""

    NO_FUNC = "NO_FUNC"
    GENERAL = "{}"
    COMPUTE = "compute_{}"
    PLOT = "plot_{}"
