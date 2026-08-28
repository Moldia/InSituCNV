from __future__ import annotations

import scanpy as sc

from insitucnv import datasets


def test_example_dataset_path_respects_env(monkeypatch, tmp_path):
    target = tmp_path / "my.h5ad"
    monkeypatch.setenv("INSITUCNV_EXAMPLE_H5AD", str(target))
    assert datasets.example_dataset_path() == target


def test_download_example_dataset_uses_local_override(monkeypatch, tiny_adata, tmp_path):
    local = tmp_path / "local.h5ad"
    tiny_adata.write(local)
    monkeypatch.setenv("INSITUCNV_EXAMPLE_H5AD", str(local))

    path = datasets.download_example_dataset()
    assert path == local
    sc.read_h5ad(path)


def test_download_example_dataset_missing_local(monkeypatch, tmp_path):
    monkeypatch.setenv("INSITUCNV_EXAMPLE_H5AD", str(tmp_path / "nope.h5ad"))
    try:
        datasets.download_example_dataset()
    except FileNotFoundError:
        pass
    else:
        raise AssertionError("expected FileNotFoundError")
