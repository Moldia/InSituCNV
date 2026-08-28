"""Helpers to fetch a small example dataset for trying InSituCNV.

The example ``.h5ad`` is not shipped in the package; it is downloaded once from a
GitHub release asset and cached locally. Override the source with the
``INSITUCNV_EXAMPLE_URL`` / ``INSITUCNV_EXAMPLE_SHA256`` environment variables, or
point ``INSITUCNV_EXAMPLE_H5AD`` at a local file to skip downloading entirely.
"""

from __future__ import annotations

import hashlib
import os
import urllib.request
from pathlib import Path

__all__ = ["download_example_dataset", "example_dataset_path"]

# Uploaded as a release asset; see scripts/make_example_dataset.py.
_DEFAULT_URL = "https://github.com/Moldia/InSituCNV/releases/download/v0.1.1/insitucnv_example.h5ad"
_DEFAULT_SHA256 = ""  # fill in when the asset is uploaded; empty disables the check
_FILENAME = "insitucnv_example.h5ad"


def _cache_dir() -> Path:
    env = os.environ.get("INSITUCNV_CACHE_DIR")
    if env:
        base = Path(env)
    else:
        xdg = os.environ.get("XDG_CACHE_HOME")
        base = Path(xdg) / "insitucnv" if xdg else Path.home() / ".cache" / "insitucnv"
    base.mkdir(parents=True, exist_ok=True)
    return base


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def example_dataset_path() -> Path:
    """Return the local path the example dataset is (or would be) cached at."""
    local = os.environ.get("INSITUCNV_EXAMPLE_H5AD")
    if local:
        return Path(local).expanduser()
    return _cache_dir() / _FILENAME


def download_example_dataset(force: bool = False) -> Path:
    """Download (once) and return the path to the example ``.h5ad``.

    Parameters
    ----------
    force
        Re-download even if the cached file already exists.
    """
    local = os.environ.get("INSITUCNV_EXAMPLE_H5AD")
    if local:
        path = Path(local).expanduser()
        if not path.exists():
            raise FileNotFoundError(f"INSITUCNV_EXAMPLE_H5AD points at a missing file: {path}")
        return path

    url = os.environ.get("INSITUCNV_EXAMPLE_URL", _DEFAULT_URL)
    expected = os.environ.get("INSITUCNV_EXAMPLE_SHA256", _DEFAULT_SHA256)
    path = example_dataset_path()

    if path.exists() and not force:
        if not expected or _sha256(path) == expected:
            return path

    tmp = path.with_suffix(path.suffix + ".part")
    try:
        urllib.request.urlretrieve(url, tmp)  # noqa: S310 - fixed https release asset
    except OSError as exc:  # pragma: no cover - network dependent
        raise RuntimeError(
            f"Could not download the example dataset from {url}. Set INSITUCNV_EXAMPLE_URL "
            "or INSITUCNV_EXAMPLE_H5AD to use a different source."
        ) from exc

    if expected:
        actual = _sha256(tmp)
        if actual != expected:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(f"Checksum mismatch for {url}: expected {expected}, got {actual}.")

    tmp.replace(path)
    return path
