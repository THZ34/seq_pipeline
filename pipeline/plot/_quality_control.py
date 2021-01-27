import pandas as pd
import seaborn as sns



def pct_counts(adata, kind='violin'):
    """"""

    if 'qc' not in adata.uns:
        import scanpy as sc
        adata.uns['qc'] = sc.pp.calculate_qc_metrics(adata)
    else:
        if 'mt_pct' not in adata.uns['qc'][0].columns:
            mt_gene = [i for i in adata.var_names if i.startswith('MT-')]
            adata.uns['qc'][0]['pct_mt'] = adata.to_df()[mt_gene].sum(1) / adata.to_df().sum(1)


    plot_df = pd.DataFrame()
    tmp_df = adata.uns['qc'][0]
    status_dict = {'pct_counts_in_top_100_genes': 'top 100 genes',
                   'pct_counts_in_top_200_genes': 'top 200 genes',
                   'pct_counts_in_top_500_genes': 'top 500 genes',
                   'pct_mt': 'mitochondria gene'}
    for col in ['pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes',
                'pct_mt']:
        tmp_df['status'] = status_dict[col]
        tmp_df_child = tmp_df[[col, 'status']]
        tmp_df_child.columns = ['Counts Percent', 'Gene set']
        plot_df = pd.concat([plot_df, tmp_df_child], axis=0)
    if kind == 'bar':
        sns.barplot(data=plot_df, x='Gene set', y='Counts Percent')
    elif kind=='violin':
        sns.violinplot(data=plot_df, x='Gene set', y='Counts Percent')
        sns.stripplot(data=plot_df, x='Gene set', y='Counts Percent')
    else:
        sns.boxplot(data=plot_df, x='Gene set', y='Counts Percent')
        sns.stripplot(data=plot_df, x='Gene set', y='Counts Percent')

