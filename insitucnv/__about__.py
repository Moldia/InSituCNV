"""Single source of truth for the package version.

``pyproject.toml`` reads ``__version__`` from here via
``[tool.setuptools.dynamic]``. Bump this string when releasing and tag the
matching ``vX.Y.Z`` on GitHub.
"""

from __future__ import annotations

__version__ = "0.1.1"
