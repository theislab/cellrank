---
substitutions:
  CI: |-
    ```{image} https://img.shields.io/github/actions/workflow/status/theislab/cellrank/test.yml?branch=main
    :alt: CI
    :target: https://github.com/theislab/cellrank/actions
    ```
  Codecov: |-
    ```{image} https://codecov.io/gh/theislab/cellrank/branch/main/graph/badge.svg
    :alt: Coverage
    :target: https://codecov.io/gh/theislab/cellrank
    ```
  Discourse: |-
    ```{image} https://img.shields.io/discourse/posts?color=yellow&logo=discourse&server=https%3A%2F%2Fdiscourse.scverse.org
    :alt: Discourse
    :target: https://discourse.scverse.org/c/ecosystem/cellrank/
    ```
  Docs: |-
    ```{image} https://img.shields.io/readthedocs/cellrank
    :alt: Documentation
    :target: https://cellrank.readthedocs.io/
    ```
  Downloads: |-
    ```{image} https://static.pepy.tech/badge/cellrank
    :alt: Downloads
    :target: https://pepy.tech/project/cellrank
    ```
  PyPI: |-
    ```{image} https://img.shields.io/pypi/v/cellrank.svg
    :alt: PyPI
    :target: https://pypi.org/project/cellrank
    ```
---

{{ PyPI }} {{ Downloads }} {{ CI }} {{ Docs }} {{ Codecov }} {{ Discourse }}

# CellRank 2: Unified fate mapping in multiview single-cell data

```{image} _static/img/light_mode_overview.png
:align: center
:class: only-light
:width: 600px
```

```{image} _static/img/dark_mode_overview.png
:align: center
:class: only-dark
:width: 600px
```

**CellRank** {cite}`lange:22,weiler:24` is a modular framework to study cellular dynamics based on Markov state modeling of
multi-view single-cell data. See {doc}`about CellRank <about/index>` to learn more and {doc}`our citation guide <about/cite>` for guidance on
citing our work correctly. Two peer-reviewed publications accompany our software:

- [CellRank for directed single-cell fate mapping](https://doi.org/10.1038/s41592-021-01346-6)
- [CellRank 2: unified fate mapping in multiview single-cell data](https://doi.org/10.1038/s41592-024-02303-9)

:::{important}
Please refer to {doc}`our citation guide <about/cite>` to cite our software correctly.
:::

CellRank scales to large cell numbers, is fully compatible with the [scverse] ecosystem, and is easy to use. In the
backend, it is powered by the [pyGPCCA package](https://github.com/msmdev/pyGPCCA) {cite}`reuter:19,reuter:22`. Feel
free to open an [issue] if you encounter a bug, need our help or just want to make a comment/suggestion.

## CellRank's key applications

- Estimate differentiation direction based on a varied number of biological priors, including
  {doc}`pseudotime <notebooks/tutorials/kernels/300_pseudotime>`,
  {doc}`developmental potential <notebooks/tutorials/kernels/400_cytotrace>`,
  {doc}`RNA velocity <notebooks/tutorials/kernels/200_rna_velocity>`,
  {doc}`experimental time points <notebooks/tutorials/kernels/500_real_time>`, and {mod}`more <cellrank.kernels>`.
- Compute initial, terminal and intermediate {doc}`macrostates <notebooks/tutorials/estimators/600_initial_terminal>`
  {cite}`reuter:19,reuter:22`.
- Infer {doc}`fate probabilities and driver genes <notebooks/tutorials/estimators/700_fate_probabilities>`.
- Visualize and cluster {doc}`gene expression trends <notebooks/tutorials/estimators/800_gene_trends>`.
- ... and much more, check out our {doc}`API <api/index>`.

## Getting started with CellRank

We have {doc}`notebooks/tutorials/index` to help you getting started. To see CellRank in action, explore our
manuscripts {cite}`lange:22,weiler:24` in Nature Methods.

## Contributing

We actively encourage any contribution! To get started, please check out the {doc}`contributing`.

```{toctree}
:caption: General
:hidden: true
:maxdepth: 3

installation
api/index
notebooks/tutorials/index
release_notes
contributing
references
```

```{toctree}
:caption: About
:hidden: true
:maxdepth: 3

about/index
about/version2
about/cite
about/team
GitHub <https://github.com/theislab/cellrank>
Discourse <https://discourse.scverse.org/c/ecosytem/cellrank/40>
```

[issue]: https://github.com/theislab/cellrank/issues/new/choose
[scverse]: https://scverse.org/
