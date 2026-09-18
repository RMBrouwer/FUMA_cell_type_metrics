# Script to further clean up h5ad files; 
# Exclude non-protein coding genes
# Remove low-expressed genes (and housekeeping genes?!) by running ANOVA's and removing all of those that to not differ amongst the celltypes

# Following approaches used in Skene/Bryois but also CELLECT


import anndata
import pandas as pd
import argparse
import scanpy as sc
import numpy as np
import scipy.sparse
from scipy import stats

def main(args):
    # read in the scRNAseq data
    adata_original = anndata.read_h5ad(args.h5ad)
    ct_column = "cell_type_" + args.level
    keep_column = "keep_" + args.level
    outname = args.output
    
    # because the scRNAseq data contains genes that don't have ensemble, it might be cleaner to create a filtered adata where:
    # genes that were not convertible to ensembl
    # adata_original = adata_original[adata_original.obs[keep_column] == True]
    adata_original.var_names_make_unique()
    gene_list = list(adata_original.var[adata_original.var["ensembl"].notnull()]['ensembl']) #this obtains a list of gene in symbol that has an ensemble conversion
    protein_coding_genes = list(pd.read_table("/projects/0/vusr0455/sumstats/magma_gene/gene_locations/magma_protein_coding_ensemble_genelist.txt",header=None)[0])
    gene_list = [x for x in gene_list if x in protein_coding_genes ]
    inx = [x for x in adata_original.var.index if adata_original.var.loc[x]["ensembl"] in gene_list]
    adata = adata_original[:, inx].copy() #subset based on the genes (keep the genes if it was converted to ensemble successfully)
    
    # remove cells belonging to small clusters
    cluster_sizes = adata.obs[ct_column].value_counts()
    min_cluster_size = 20
    small_clusters = cluster_sizes[cluster_sizes < min_cluster_size].index
    cells_to_keep = ~adata.obs[ct_column].isin(small_clusters)
    adata = adata[cells_to_keep].copy()

    # Run anova's using the cell type label as grouping factor on normalised data (no log-transform in the Skene paper).
    sc.pp.normalize_total(adata, target_sum=1e6)
    # borrowed this code from cellex https://github.com/perslab/CELLEX/blob/master/cellex/preprocessing/anova.py

    idx = [np.where(adata.obs[ct_column] == i)[0] for i in pd.unique(adata.obs[ct_column])]

    ### Iterate over genes and compute ANOVA for each gene
    def _anova_(row, groups):
        """Apply ANOVA to a vector
        Split the annotation by
        """
        fval, pval = stats.f_oneway(*[row[g] for g in groups])
        return np.array([fval, pval])

    ### Create dataframe for ANOVA results
    f = np.apply_along_axis(func1d=_anova_, axis=0, arr=adata.to_df(), groups=(idx))

    df_anova = pd.DataFrame({"statistic": f[0,:],"pvalue": f[1,:]},index=adata.var.index)
    idx_thresh = (df_anova.values[:,1] < 0.00001 )
    genelist_final = adata.var['ensembl'][idx_thresh]
    adata_final = adata_original[:,genelist_final.index].copy()
    
    adata_final.write_h5ad(outname)

def parse_args():
    parser = argparse.ArgumentParser(description='Generate input for extra filtering of scRNAseq data')
    parser.add_argument('--h5ad', required=True, help='Path to the h5ad data. This is the file after the "preprocessing" step.')
    parser.add_argument('--level', required=True, help='Cell type resolution level.')
    parser.add_argument('--output', required=True, help='Name of output file.')
    return parser.parse_args()



main(parse_args())
