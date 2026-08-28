Quickstart
==========

Prepare an ``AnnData`` object with raw counts (in ``layers["raw_counts"]`` or
``X``), spatial coordinates in ``obsm["spatial"]``, human gene symbols in
``var_names``, and an ``obs`` column identifying non-tumor reference cells. A
neighbor graph and genomic gene coordinates are added automatically.

.. code-block:: python

   import scanpy as sc
   import insitucnv as icv

   # Try it with the bundled example, or use sc.read_h5ad("sample.h5ad")
   adata = sc.read_h5ad(icv.download_example_dataset())

   result = icv.run_insitucnv(
       adata,
       output_dir="results/sample",
       reference_key="cell_type",
       reference_categories=["T_cells", "B_cells", "Stromal"],
       raw_layer="raw_counts",
       cluster_resolutions=[0.1, 0.2, 0.3],
       # build_neighbors=True by default; gene_reference_path=... for a custom
       # gene_name,chromosome,start,end table
   )

   cnv_adata = result["adata"]
   primary_cluster_key = result["primary_cluster_key"]

Review the generated heatmaps and spatial plots before manually assigning
normal, tumor, or clone labels with :func:`insitucnv.tl.assign_cnv_status`.
