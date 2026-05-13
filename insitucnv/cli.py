from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

from insitucnv.pipeline import annotate_cell_types, load_xenium_dataset, preprocess_expression
from insitucnv.workflow import run_insitucnv


def _float_list(value: str) -> list[float]:
    return [float(item) for item in value.split(",") if item.strip()]


def _str_list(value: str | None) -> list[str] | None:
    if value is None or value.strip() == "":
        return None
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description="Run InSituCNV on a prepared AnnData or Xenium output directory.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--input-dir", help="Xenium output directory with cell_feature_matrix.h5 and cells.csv.gz.")
    group.add_argument("--input-h5ad", help="Prepared AnnData file.")
    parser.add_argument("--sample-id", help="Sample name to write into adata.obs['sample'].")
    parser.add_argument("--output-dir", required=True, help="Directory where outputs will be written.")
    parser.add_argument("--annotation-csv", help="Optional CSV with external cell-type annotations.")
    parser.add_argument("--reference-key", default="cell_type", help="obs column used as inferCNV reference labels.")
    parser.add_argument("--reference-categories", help="Comma-separated reference categories. Auto-detected if omitted.")
    parser.add_argument("--raw-layer", default="raw_counts", help="Layer containing raw counts.")
    parser.add_argument("--target-sum", type=float, default=1e4, help="Target sum for normalization.")
    parser.add_argument("--min-counts", type=int, default=20, help="Minimum counts per cell during Xenium preprocessing.")
    parser.add_argument("--min-genes", type=int, default=10, help="Minimum detected genes per cell during Xenium preprocessing.")
    parser.add_argument("--neighbors", type=int, default=15, help="Neighbors for expression graph construction.")
    parser.add_argument("--smooth-neighbors", type=int, default=100, help="Neighbors used for smoothing before inferCNV.")
    parser.add_argument("--window-size", type=int, default=60, help="infercnvpy window size.")
    parser.add_argument("--step", type=int, default=10, help="infercnvpy step size.")
    parser.add_argument("--lfc-clip", type=float, default=4.0, help="infercnvpy lfc_clip value.")
    parser.add_argument("--chunksize", type=int, default=1000, help="infercnvpy chunksize.")
    parser.add_argument("--cluster-resolutions", default="0.1,0.2,0.3", help="Comma-separated Leiden resolutions.")
    parser.add_argument("--primary-resolution", type=float, help="Resolution used for cnv_status annotation.")
    parser.add_argument("--select-resolution-by-metrics", action="store_true", help="Select primary resolution by quality metrics.")
    parser.add_argument("--evaluate-resolution-metrics", action="store_true", help="Save quality metrics for all resolutions.")
    parser.add_argument("--run-umap", action="store_true", help="Compute and save a CNV UMAP plot.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)
    sample_id = args.sample_id

    if args.input_dir:
        adata = load_xenium_dataset(args.input_dir, sample_id=sample_id)
        adata = preprocess_expression(
            adata,
            min_counts=args.min_counts,
            min_genes=args.min_genes,
            n_neighbors=args.neighbors,
        )
        adata = annotate_cell_types(
            adata,
            annotation_csv=args.annotation_csv,
            cluster_key="leiden",
            output_key=args.reference_key,
        )
    else:
        adata = sc.read_h5ad(args.input_h5ad)
        if sample_id is not None:
            adata.obs["sample"] = sample_id
        if "cell_id" not in adata.obs:
            adata.obs["cell_id"] = adata.obs_names.astype(str)
        if args.raw_layer not in adata.layers:
            adata.layers[args.raw_layer] = adata.X.copy()

    run_insitucnv(
        adata,
        output_dir=Path(args.output_dir),
        reference_key=args.reference_key,
        reference_categories=_str_list(args.reference_categories),
        raw_layer=args.raw_layer,
        target_sum=args.target_sum,
        smoothing_neighbors=args.smooth_neighbors,
        window_size=args.window_size,
        step=args.step,
        lfc_clip=args.lfc_clip,
        chunksize=args.chunksize,
        cluster_resolutions=_float_list(args.cluster_resolutions),
        primary_resolution=args.primary_resolution,
        select_resolution_by_metrics=args.select_resolution_by_metrics,
        evaluate_resolution_metrics=args.evaluate_resolution_metrics,
        run_umap=args.run_umap,
    )


if __name__ == "__main__":
    main()
