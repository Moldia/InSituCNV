from __future__ import annotations

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

from insitucnv import pl, pp, tl

try:
    __version__ = version("insitucnv")
except PackageNotFoundError:  # pragma: no cover - local tree without installation
    __version__ = "0.0.0"

_PIPELINE_EXPORTS = {
    "annotate_cell_types",
    "load_xenium_dataset",
    "preprocess_expression",
    "run_xenium_cnv_protocol",
}
_WORKFLOW_EXPORTS = {"run_insitucnv"}


def __getattr__(name: str):
    if name in _PIPELINE_EXPORTS:
        module = import_module("insitucnv.pipeline")
    elif name in _WORKFLOW_EXPORTS:
        module = import_module("insitucnv.workflow")
    else:
        raise AttributeError(f"module 'insitucnv' has no attribute {name!r}")

    value = getattr(module, name)
    globals()[name] = value
    return value

__all__ = [
    "__version__",
    "annotate_cell_types",
    "load_xenium_dataset",
    "pl",
    "pp",
    "preprocess_expression",
    "run_insitucnv",
    "run_xenium_cnv_protocol",
    "tl",
]
