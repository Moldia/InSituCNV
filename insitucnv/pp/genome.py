from __future__ import annotations

from pathlib import Path

import pandas as pd

GENOMIC_COLUMNS = ("chromosome", "start", "end")


def _read_reference_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _default_reference() -> pd.DataFrame:
    """Return the infercnvpy reference gene table used in the manuscript notebooks."""
    import infercnvpy as cnv

    return cnv.datasets.maynard2020_3k().var.copy()


def add_genomic_positions(
    adata,
    reference: pd.DataFrame | None = None,
    reference_path: str | Path | None = None,
    gene_name_col: str = "gene_name",
    drop_unmapped_genes: bool = True,
    copy: bool = True,
    quiet: bool = False,
):
    """Add genomic position columns required by ``infercnvpy``.

    Parameters
    ----------
    adata
        AnnData object whose ``var_names`` are gene symbols.
    reference
        Optional gene annotation table. It must contain ``gene_name_col`` plus
        ``chromosome``, ``start``, and ``end`` columns. If omitted, the
        infercnvpy Maynard 2020 gene reference is used, matching the manuscript
        notebooks.
    reference_path
        Optional CSV/TSV gene annotation table. Ignored when ``reference`` is
        provided.
    gene_name_col
        Column in the reference table matching ``adata.var_names``.
    drop_unmapped_genes
        Drop genes that have no genomic coordinates. ``infercnvpy`` requires
        genomic coordinates, so this should usually stay ``True``.
    copy
        Return a copy. When ``drop_unmapped_genes=True`` this function returns a
        subset AnnData object, so assigning the result is recommended.
    quiet
        Suppress the short mapping summary.

    Returns
    -------
    AnnData
        AnnData with ``chromosome``, ``start``, and ``end`` in ``adata.var``.
    """
    if reference is None:
        reference = _read_reference_table(reference_path) if reference_path is not None else _default_reference()
    else:
        reference = reference.copy()

    missing = [col for col in (gene_name_col, *GENOMIC_COLUMNS) if col not in reference.columns]
    if missing:
        raise KeyError(f"Reference table is missing required column(s): {missing}")

    reference[gene_name_col] = reference[gene_name_col].astype(str)
    genes_total = adata.n_vars
    mapped = adata.var_names.astype(str).isin(reference[gene_name_col])

    if drop_unmapped_genes:
        out = adata[:, mapped].copy()
    else:
        out = adata.copy() if copy else adata

    reference = reference.drop_duplicates(gene_name_col).set_index(gene_name_col)
    for col in ("gene_id", "symbol", "ensg", "chromosome", "start", "end"):
        if col in reference.columns:
            out.var[col] = out.var_names.map(reference[col])

    if "symbol" not in out.var.columns:
        out.var["symbol"] = out.var_names

    genes_located = int(out.var["chromosome"].notna().sum())
    if not quiet:
        print(f"Added genomic positions to {genes_located} genes out of {genes_total}.")

    return out
