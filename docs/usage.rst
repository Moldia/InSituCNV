Usage
=====

InSituCNV is designed to be used both as a Python package and through provided notebooks.

Python API
----------

You can import InSituCNV in your Python scripts or notebooks:

.. code-block:: python

   import insitucnv as icv

   # Example: smooth data for CNV inference
   icv.tl.smooth_data_for_cnv(adata)

   # Example: run inferCNV
   icv.tl.run_infercnv(adata, reference_key='cell_type', reference_categories=['Stromal'])

For more detailed information, see the Tutorials section.
