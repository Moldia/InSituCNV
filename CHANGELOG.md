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
  requiring a precomputed graph (which must now be built with
  `smoothing_neighbors`, see Fixed).
- `add_genomic_positions()` keeps `chromosome/start/end` already present in
  `adata.var` instead of always consulting the (network) reference table.
- `scripts/synthetic_example_dataset.py` builds a small offline dataset used by
  the notebook-execution CI job.
- `CONTRIBUTING.md`, this changelog, and an (unenforced) ruff config.

### Changed
- Single high-level engine: `run_xenium_cnv_protocol()` now prepares an AnnData
  and delegates to `run_insitucnv()` instead of re-implementing the smoothing /
  inferCNV / clustering steps. It returns `(adata, summary)` (0.1.0 returned
  `(adata, metrics, summary)`).
- `smoothing_neighbors` defaults to `20` everywhere (`smooth_data_for_cnv()`,
  `prepare_cnv_input()`, `run_insitucnv()` and the CLI).
- The shipped notebook `notebooks/run_insitucnv.ipynb` was rewritten to use the
  package API end to end (no direct `infercnvpy` calls) and to run top to bottom
  on the example dataset.
- Package version is now single-sourced from `insitucnv/__about__.py`.

### Removed
- `neighbors_n_neighbors` of `prepare_cnv_input()` / `run_insitucnv()`: the
  smoothing graph is built with `smoothing_neighbors`. The CLI `--neighbors` now
  only sets the expression graph used for cell-type clustering.
- `insitucnv.analysis.find_optimal_clustering()` and the
  `select_resolution_by_metrics` / `evaluate_resolution_metrics` options of
  `run_insitucnv()` (and the matching CLI flags). Pick a clustering resolution by
  inspecting the chromosome heatmap and pass `primary_resolution=`.

### Fixed
- `smoothing_neighbors` now actually sets the smoothing neighborhood. Previously,
  when it was larger than the `n_neighbors` of the stored graph, scVelo only
  warned and smoothing silently used the whole smaller graph (e.g. the default
  15-neighbor graph when 20 or 100 were requested). Now `prepare_cnv_input()` /
  `run_insitucnv()` reuse `adata.uns["neighbors"]` only if it was built with
  exactly `smoothing_neighbors`, and otherwise build a dedicated graph under
  `adata.uns["cnv_smoothing"]` (no longer overwriting the user's graph or
  `X_pca`); `smooth_data_for_cnv()` gains `neighbors_key` and raises if the graph's
  `n_neighbors` differs from the requested one. **Results with default settings
  change**, and earlier runs are not reproduced by passing the same
  `smoothing_neighbors`.
- `cluster_cnv_resolutions()` and the Xenium spatial plots no longer use the
  `matplotlib.cm.get_cmap` / `plt.get_cmap(name, N)` APIs removed in
  matplotlib 3.9, so clustering works on modern matplotlib.
- `cluster_cnv_resolutions()` and `run_insitucnv()` share one CNV Leiden key
  convention (`cnv_leiden_res<r>`) via `insitucnv.tl.cnv_leiden_key()`.

## [0.1.0] - 2026-05-23

- Initial packaged release of the InSituCNV workflow from Jensen et al.
  (bioRxiv 2025, doi:10.1101/2025.07.02.662761).
