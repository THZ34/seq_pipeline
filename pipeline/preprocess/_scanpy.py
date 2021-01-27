import scanpy as sc


def sc_recipe(adata):
    #
    sc.pl.highest_expr_genes(adata, n_top=20)
    adata.obs['n_counts'] = adata.to_df().sum(1)
    adata.obs['n_genes'] = adata.to_df().astype(dtype=bool).astype(dtype=int).sum(1)
    mito_gene = [i for i in adata.var_names if i.startswith('MT-')]
    adata.obs['mito-percent'] = adata.to_df()[mito_gene].sum(1) / adata.to_df().sum(1)


