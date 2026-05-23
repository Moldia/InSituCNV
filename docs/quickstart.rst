Quickstart
==========

Prepare an ``AnnData`` object with raw counts, a neighbor graph, genomic gene
coordinates, and an ``obs`` column identifying normal reference cells.

.. code-block:: python

   import scanpy as sc
   import insitucnv as icv

   adata = sc.read_h5ad("sample.h5ad")

   result = icv.run_insitucnv(
       adata,
       output_dir="results/sample",
       reference_key="cell_type",
       reference_categories=["T_cells", "B_cells", "Stromal"],
       raw_layer="raw_counts",
       cluster_resolutions=[0.1, 0.2, 0.3],
   )

   cnv_adata = result["adata"]
   primary_cluster_key = result["primary_cluster_key"]

Review the generated heatmaps and spatial plots before manually assigning
normal, tumor, or clone labels.
