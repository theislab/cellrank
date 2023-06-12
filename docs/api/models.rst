.. module:: cellrank.models

Models
======
Models fit gene expression trends in pseudotime; they assume some parametric form for the gene trend and estimate
parameters using an objective function. Note that some models require you to have R and rpy2 installed.

.. currentmodule:: cellrank

.. autosummary::
    :toctree: _autosummary/models

    models.GAM
    models.GAMR
    models.SKLearnModel
