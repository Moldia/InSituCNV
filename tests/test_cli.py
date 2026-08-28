from __future__ import annotations

from insitucnv.cli import main


def test_cli_input_h5ad(tiny_adata, gene_reference_csv, tmp_path):
    h5ad = tmp_path / "input.h5ad"
    tiny_adata.write(h5ad)
    out_dir = tmp_path / "out"

    main(
        [
            "--input-h5ad",
            str(h5ad),
            "--output-dir",
            str(out_dir),
            "--reference-key",
            "cell_type",
            "--reference-categories",
            "T_cells",
            "--gene-reference-path",
            gene_reference_csv,
            "--cluster-resolutions",
            "0.5,1.0",
            "--window-size",
            "10",
            "--step",
            "5",
            "--smooth-neighbors",
            "15",
        ]
    )

    assert (out_dir / "adata_cnv.h5ad").exists()
    assert (out_dir / "run_summary.json").exists()
