# Contributing to InSituCNV

Thanks for your interest in improving InSituCNV. This document covers the
development workflow. `pyproject.toml` is the single source of truth for runtime
dependencies; `insitucnv.yml` is a convenience conda environment.

## Development setup

```bash
git clone https://github.com/Moldia/InSituCNV.git
cd InSituCNV
conda env create -f insitucnv.yml
conda activate insitucnv_env
pip install -e ".[dev,docs]"
pre-commit install
```

Notebooks and examples import the installed package directly (`import insitucnv as
icv`); no `sys.path` edits are needed with the editable install.

## Tests

```bash
pytest -m "not network" --cov=insitucnv
```

`network`-marked tests download reference data and are skipped by default. Run the
notebook end to end with:

```bash
pytest --nbmake notebooks/run_insitucnv.ipynb
```

## Linting and formatting

```bash
ruff check .
ruff format --check .
```

`pre-commit` runs these automatically on commit.

## Documentation

```bash
sphinx-build -b html docs docs/_build/html
```

## Releasing

1. Move the `Unreleased` section of `CHANGELOG.md` under a new `## [X.Y.Z]` heading.
2. Bump `__version__` in `insitucnv/__about__.py` to match.
3. Commit, then tag: `git tag vX.Y.Z && git push --tags`.
4. Create a GitHub Release for the tag. Publishing to PyPI happens automatically
   via Trusted Publishing (`.github/workflows/publish-pypi.yml`).

Do not store PyPI or TestPyPI API tokens in the repository. TestPyPI releases can
be triggered manually from the `publish-testpypi.yml` workflow.
