Installation
============

InSituCNV requires Python 3.10 or later. Install the released package from
PyPI when it is available:

.. code-block:: bash

   pip install insitucnv

For development, install the package from a local checkout in editable mode:

.. code-block:: bash

   git clone https://github.com/Moldia/InSituCNV.git
   cd InSituCNV
   pip install -e ".[dev,docs]"

To build the documentation locally:

.. code-block:: bash

   sphinx-build -b html docs docs/_build/html

To build source and wheel distributions:

.. code-block:: bash

   python -m build

Requirements
------------

InSituCNV requires Python 3.10 or later and depends on several key packages:

* `scanpy <https://scanpy.readthedocs.io/>`_
* `infercnvpy <https://infercnvpy.readthedocs.io/>`_
* `scvelo <https://scvelo.readthedocs.io/>`_
* `anndata <https://anndata.readthedocs.io/>`_
