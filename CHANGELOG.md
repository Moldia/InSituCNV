# Changelog

All notable changes to InSituCNV are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `insitucnv.download_example_dataset()` fetches a small demo `.h5ad` so the
  notebook and quickstart run out of the box without preparing your own data.
- `prepare_cnv_input()` / `run_insitucnv()` gain `build_neighbors` (default
  `True`): when the AnnData has no neighbor graph, the standard
  normalize → log1p → PCA → `sc.pp.neighbors` graph used for smoothing is built
  automatically. Pass `build_neighbors=False` for the previous behaviour of
  requiring a precomputed graph.
- `add_genomic_positions()` keeps `chromosome/start/end` already present in
  `adata.var` instead of always consulting the (network) reference table.
- `scripts/synthetic_example_dataset.py` builds a small offline dataset used by
  the notebook-execution CI job.
- `CONTRIBUTING.md`, this changelog, a ruff config and a pre-commit config.

### Changed
- Single high-level engine: `run_xenium_cnv_protocol()` now prepares an AnnData
  and delegates to `run_insitucnv()` instead of re-implementing the smoothing /
  inferCNV / clustering steps. Its `smoothing_neighbors` default changed from
  `20` to `100` to match `run_insitucnv()`. **This can change results for the
  Xenium protocol path** — pass `smoothing_neighbors=20` to reproduce older runs.
- `find_optimal_clustering()` now writes cluster keys as `cnv_leiden_res<r>`,
  matching `cluster_cnv_resolutions()` and `select_best_resolution()` (previously
  `cnv_leiden_<r>`, which broke `select_resolution_by_metrics=True`).
- The shipped notebook `notebooks/run_insitucnv.ipynb` was rewritten to use the
  package API end to end (no direct `infercnvpy` calls) and to run top to bottom
  on the example dataset.
- Package version is now single-sourced from `insitucnv/__about__.py`.

### Fixed
- `cluster_cnv_resolutions()` and the Xenium spatial plots no longer use the
  `matplotlib.cm.get_cmap` / `plt.get_cmap(name, N)` APIs removed in
  matplotlib 3.9, so clustering works on modern matplotlib.

## [0.1.0] - 2026-05-23

- Initial packaged release of the InSituCNV workflow from Jensen et al.
  (bioRxiv 2025, doi:10.1101/2025.07.02.662761).
