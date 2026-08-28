from datetime import datetime
from importlib.metadata import PackageNotFoundError, version as package_version

# -- Project information -----------------------------------------------------
project = "InSituCNV"
copyright = f"{datetime.now().year}, Augusta Jensen"
author = "Augusta Jensen"

try:
    release = package_version("insitucnv")
except PackageNotFoundError:
    release = "0.0.0"
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_nb",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

autosummary_generate = True
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = True

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# -- Options for intersphinx extension ---------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
    "infercnvpy": ("https://infercnvpy.readthedocs.io/en/latest/", None),
    "scvelo": ("https://scvelo.readthedocs.io/en/latest/", None),
}

# -- Notebook options --------------------------------------------------------
nb_execution_mode = "off"
nb_execution_timeout = 120
myst_heading_anchors = 3
