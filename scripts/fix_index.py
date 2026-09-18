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
import os 

def touch(fname):
    if os.path.exists(fname):
        os.utime(fname, None)
    else:
        open(fname, 'a').close()

def main(args):
    # read in the scRNAseq data
    adata_original = anndata.read_h5ad(args.h5ad)
    outname = args.output
    outname_donotrun = args.output[:-4] + 'norun.touch'
    keep_level = "keep_" + args.level
    ct_level = "cell_type_" + args.level
 
    # because the scRNAseq data contains genes that don't have ensemble, it might be cleaner to create a filtered adata where:
    # genes that were not convertible to ensembl
    # adata_original = adata_original[adata_original.obs[keep_column] == True]
    adata_original.var_names_make_unique()
    gene_list = list(adata_original.var[adata_original.var["ensembl"].notnull()]['ensembl']) #this obtains a list of gene in symbol that has an ensemble conversion
    inx = [x for x in adata_original.var.index if adata_original.var.loc[x]["ensembl"] in gene_list]
    adata = adata_original[:, inx].copy() #subset based on the genes (keep the genes if it was converted to ensemble successfully)
    
    # In these datasets, we have some issues with duplicate ensembl ids. These stem from duplicate symbols with different 'other' identifiers such as ensembl ids. 
    # We chose to re-match the symbols to ensembl ids even if these were present to have every dataset be based on the same conversion file
    # Now we need to remove duplicate ids
    # This seems a bit complicated - in principle we should trust the original gene naming more than our own; however, duplicates arise mostly for non-coding or pseudo genes and these are not used in subsequent analyses. We expect the overall effect of this to be minimal. 
    # For the Siletti data - we can make use of the fact that original Ids are present, we remove the duplicates for which the ids do not match. In any other case, we remove both duplicates
    
    if 'ensembl_ids' in adata.var.columns:
        inx_remove = [x for x in adata.var.index if adata.var.loc[x]['ensembl'] != adata.var.loc[x]['ensembl_ids'] and x in adata.var[adata.var['ensembl'].duplicated(keep=False)].index]
        inx = [x for x in adata.var.index if x not in inx_remove] 
        adata = adata[:,inx].copy()
    # remove any remaining duplicates
    inx_remove = [x for x in adata.var[adata.var['ensembl'].duplicated(keep=False)].index]
    inx = [x for x in adata.var.index if x not in inx_remove]
    adata = adata[:,inx].copy()

    # In case there is a keep_level variable present, remove the cells at that level
    if keep_level in adata.obs.columns:
        adata = adata[adata.obs[keep_level] == True]
    adata.var.index = adata.var['ensembl'].astype(str)
    # remove the cells without label at this level
    adata = adata[~adata.obs[ct_level].isna()]
    if adata.obs[ct_level].value_counts().max() < 20:
        touch(outname_donotrun)
    else: 
        adata.write_h5ad(outname)

def parse_args():
    parser = argparse.ArgumentParser(description='Generate filtered input in which missing ensembl names are removed and the index = ensembl')
    parser.add_argument('--h5ad', required=True, help='Path to the h5ad data. This is the file after the "preprocessing" step.')
    parser.add_argument('--output', required=True, help='Name of output file.')
    parser.add_argument('--level', required=True, help='Level to filter cells on')
    return parser.parse_args()



main(parse_args())
