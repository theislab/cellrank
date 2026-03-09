[![PyPI](https://img.shields.io/pypi/v/cellrank.svg)](https://pypi.org/project/cellrank)
[![Downloads](https://static.pepy.tech/badge/cellrank)](https://pepy.tech/project/cellrank)
[![CI](https://img.shields.io/github/actions/workflow/status/theislab/cellrank/test.yaml?branch=main)](https://github.com/theislab/cellrank/actions)
[![Docs](https://img.shields.io/readthedocs/cellrank)](https://cellrank.readthedocs.io/)
[![Codecov](https://codecov.io/gh/theislab/cellrank/branch/main/graph/badge.svg)](https://codecov.io/gh/theislab/cellrank)
[![Discourse](https://img.shields.io/discourse/posts?color=yellow&logo=discourse&server=https%3A%2F%2Fdiscourse.scverse.org)](https://discourse.scverse.org/c/ecosystem/cellrank/)
[![Python Version](https://img.shields.io/pypi/pyversions/cellrank)](https://pypi.org/project/cellrank)

# CellRank 2: Unified fate mapping in multiview single-cell data

<!-- image-start -->
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/scverse/cellrank/main/docs/_static/img/dark_mode_overview.png">
  <img alt="CellRank overview" src="https://raw.githubusercontent.com/scverse/cellrank/main/docs/_static/img/light_mode_overview.png" width="600px" align="center">
</picture>
<!-- image-end -->

**CellRank** is a modular framework to study cellular dynamics based on Markov state modeling of
multi-view single-cell data. See our [documentation], and the [CellRank 1] and [CellRank 2] manuscripts to learn more.
Read a summary of the CellRank papers [here](https://cellrank.readthedocs.io/en/latest/about/cite.html#cellrank-papers).

⚠️ **Please refer to [our citation guide](https://cellrank.readthedocs.io/en/latest/about/cite.html) to cite our software correctly.**

CellRank scales to large cell numbers, is fully compatible with the [scverse] ecosystem, and is easy to use.
In the backend, it is powered by [pyGPCCA] ([Reuter et al. (2018)]). Feel
free to open an [issue] if you encounter a bug, need our help, or just want to make a comment/suggestion.

## CellRank's key applications

- Estimate differentiation direction based on a varied number of biological priors, including RNA velocity
  ([La Manno et al. (2018)], [Bergen et al. (2020)]), any pseudotime or developmental potential,
  experimental time points, metabolic labels, and more.
- Compute initial, terminal and intermediate macrostates.
- Infer fate probabilities and driver genes.
- Visualize and cluster gene expression trends.
- ... and much more, check out our [documentation].

## Installation

```bash
pip install cellrank
```

See the [installation guide](https://cellrank.readthedocs.io/en/latest/installation.html) for more options.

## Related packages

If you like CellRank, check out these packages from the same authors.
Almost all are part of the [scverse ecosystem].

| Package | Description | Reference |
|---------|-------------|-----------|
| [moscot] | Optimal transport for temporal, spatial, and spatio-temporal single-cell mapping | [Klein et al. (2025)] |
| [moslin] | Trajectory inference with lineage barcodes via optimal transport (part of moscot) | [Lange et al. (2024)] |
| [VeloVI] | RNA velocity with variational inference and uncertainty quantification (part of scvi-tools) | [Gayoso et al. (2024)] |
| [RegVelo] | Jointly learning gene regulation and RNA velocity | [Wang et al. (2024)] |
| [CellMapper] | kNN-based label, embedding, and molecular layer transfer between datasets | — |
| [CellAnnotator] | LLM-based cell type annotation with support for major LLM providers | — |

[moscot]: https://moscot.readthedocs.io/
[moslin]: https://moscot.readthedocs.io/en/latest/notebooks/tutorials/100_lineage.html
[VeloVI]: https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/velovi.html
[RegVelo]: https://regvelo.readthedocs.io/
[CellMapper]: https://cellmapper.readthedocs.io/
[CellAnnotator]: https://cell-annotator.readthedocs.io/
[scverse ecosystem]: https://scverse.org/packages/#ecosystem
[Klein et al. (2025)]: https://doi.org/10.1038/s41586-024-08453-2
[Lange et al. (2024)]: https://doi.org/10.1186/s13059-024-03422-4
[Gayoso et al. (2024)]: https://doi.org/10.1038/s41592-023-01994-w
[Wang et al. (2024)]: https://doi.org/10.1101/2024.12.11.627935

[La Manno et al. (2018)]: https://doi.org/10.1038/s41586-018-0414-6
[Bergen et al. (2020)]: https://doi.org/10.1038/s41587-020-0591-3
[Reuter et al. (2018)]: https://doi.org/10.1021/acs.jctc.8b00079
[scverse]: https://scverse.org/
[pyGPCCA]: https://github.com/msmdev/pyGPCCA
[CellRank 1]: https://www.nature.com/articles/s41592-021-01346-6
[CellRank 2]: https://doi.org/10.1038/s41592-024-02303-9
[documentation]: https://cellrank.readthedocs.io/en/latest/
[issue]: https://github.com/theislab/cellrank/issues/new/choose
