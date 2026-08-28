from __future__ import annotations

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

from insitucnv import pl, pp, tl
from insitucnv.__about__ import __version__ as _about_version

try:
    __version__ = version("insitucnv")
except PackageNotFoundError:  # pragma: no cover - local tree without installation
    __version__ = _about_version

_PIPELINE_EXPORTS = {
    "annotate_cell_types",
    "load_xenium_dataset",
    "preprocess_expression",
    "run_xenium_cnv_protocol",
}
_WORKFLOW_EXPORTS = {"run_insitucnv"}
_DATASET_EXPORTS = {"download_example_dataset", "example_dataset_path"}


def __getattr__(name: str):
    if name in _PIPELINE_EXPORTS:
        module = import_module("insitucnv.pipeline")
    elif name in _WORKFLOW_EXPORTS:
        module = import_module("insitucnv.workflow")
    elif name in _DATASET_EXPORTS:
        module = import_module("insitucnv.datasets")
    else:
        raise AttributeError(f"module 'insitucnv' has no attribute {name!r}")

    value = getattr(module, name)
    globals()[name] = value
    return value


__all__ = [
    "__version__",
    "annotate_cell_types",
    "download_example_dataset",
    "example_dataset_path",
    "load_xenium_dataset",
    "pl",
    "pp",
    "preprocess_expression",
    "run_insitucnv",
    "run_xenium_cnv_protocol",
    "tl",
]
