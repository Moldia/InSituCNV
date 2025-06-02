import infercnvpy as cnv

def add_genomic_positions(adata):
    """
    Adds genomic positions to the gene names (adata.var_names), including chromosome, start, and end. Also adds the gene id and symbol according to the dataset: cnv.datasets.maynard2020_3k().
    
    Parameters:
        adata (AnnData object): The AnnData object where genomic positions are to be added

    Output:
        adata (AnnData object): The AnnData object where genomic positions are to be added
    """

    # Add genomic positions
    dat = cnv.datasets.maynard2020_3k()
    dat.var.loc[:, ["ensg", "chromosome", "start", "end"]].head()
    genes_total = len(adata.var_names)
    adata=adata[:,adata.var.index.isin(dat.var['gene_name'])]
    genes_located = len(adata.var_names)
    
    name2gid=dict(zip(dat.var['gene_name'],dat.var['gene_id']))
    name2symbol=dict(zip(dat.var['gene_name'],dat.var.index))

    adata.var['gene_id']=adata.var.index.map(name2gid)
    adata.var['symbol']=adata.var.index.map(name2symbol)

    name2ensg=dict(zip(dat.var['gene_name'],dat.var['ensg']))
    name2chromosome=dict(zip(dat.var['gene_name'],dat.var['chromosome']))
    name2start=dict(zip(dat.var['gene_name'],dat.var['start']))
    name2end=dict(zip(dat.var['gene_name'],dat.var['end']))

    adata.var['ensg']=adata.var.index.map(name2ensg)
    adata.var['chromosome']=adata.var.index.map(name2chromosome)
    adata.var['start']=adata.var.index.map(name2start)
    adata.var['end']=adata.var.index.map(name2end)

    print(f'Added genomic positions and ids to {genes_located} genes, out of {genes_total} genes in total!')

    return adata