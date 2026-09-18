# FUMA cell type metrics
Rachel Brouwer - September 18, 2026

This repository contains the documentation and scripts to compute four different metrics for MAGMA-based GWAS-to-cell-type analysis, as implemented in the FUMA platform available at fuma.ctglab.nl.

----

# Requirements: 
snakemake (my installation ran on version 7.25.4) and several conda environments. Json files with the specifics can be found in the directory conda_envs. 

----

# scRNA preprocessing 
We assume that the datafiles are in the following format (see also https://github.com/vu-ctg/master_scripts/tree/main/Preprocessing_scRNA)

- h5ad annotated data file
- adata.X is a matrix indexed by obs x var in which var is the gene and obs is the cell
- adata.obs contains a column describing the clustering, potentially at different levels, named "cell_type_level_1", "cell_type_level_2", etc.
- adata.obs contains a column per cell type label, that describes whether a cell should be in/excluded at that level (potentially different at different levels)
- adata.var has ENSG gene names as an index - genes without corresponding ensemble gene ID should be filtered out

We assume that each dataset has at least 3 clusters, with 20 cells each. In the pipeline, smaller clusters will be filtered out, which could lead to errors if the remaining number of clusters is too small. To check whether some of the datafiles need to be excluded, you can run:

```
snakemake -s filter_datasets.smk -j 1 --configfile configfile_metrics_FUMA.json
```

Output will be a text file containing the IDs that should be removed based on this criterion: excluded_datasets_no_clusters.txt
  
----

# Running the pipeline

Adapt the json file to contain the identifiers of the scRNA datafiles that you want to run. 
The most efficient way to run this pipeline depends on the size of the scRNA dataset in relation to the computing cluster. Not all steps in this pipeline increase in memory use in the same way. Have a look at the "multipliers" in the "resources"-parameter in each snakemake rule. 
Submit the following to the queue for many but small datasets
```
snakemake -s create_metrics.smk -j 10 --use-conda --configfile configfile_metrics_FUMA.json --rerun-trigger mtime
```
or alternatively, create a file with "dataset_ids_input.txt", one ID per line, and run sequentially so that you can use all available memory. 
```
while read line
do
sed "s/ID/${line}/g" < config_metrics_FUMA_BASE.json > config_metrics_FUMA_temp.json
snakemake -s create_metrics.smk --cores 10 --rerun-triggers mtime --resources mem_mb=224000 --use-conda -k --configfile config_metrics_FUMA_temp.json 
done < dataset_ids_input.txt
```

---- 

# Output

Files to use as input for the MAGMA gene-property analyses: 
```
fuma/${ID}/${level}/${ID}.${level}_means_cell_log_counts_pM.tsv
cellex/${ID}/${level}/${ID}.${level}.esmu_fmt.tsv
```
Files to use as input for the MAGMA gene-set analyses:
```
ewce/${ID}/${level}/${ID}.${level}_spec_ewce_top10.txt
cepo/${ID}/${level}/${ID}.${level}_cepo_top10.txt
```
----
