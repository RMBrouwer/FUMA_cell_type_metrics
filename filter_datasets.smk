import os
import anndata

levels_per_dataset=dict()
all_levels = []
remove_list = []

for i in config["scrnas"]:
   adata = anndata.read_h5ad(config["input"] + "/" + i + "/" + i + ".h5ad",backed='r')
   levels_per_dataset[i]=[col[10:] for col in adata.obs.columns if 'cell_type_level' in col]
   for j in levels_per_dataset[i]:
      keep_level = "keep_" + j
      ct_level = "cell_type_" + j
   
   obs=adata.obs.copy() 
   # In case there is a keep_level variable present, remove the cells at that level
   if keep_level in adata.obs.columns:
      obs = obs[obs[keep_level] == True]
   # remove the cells without label at this level
   obs = obs[~obs[ct_level].isna()]
   if sum(obs[ct_level].value_counts() > 20) < 3:
      remove_list.append(i + "_" + j)
   

with open('exclude_datasets_no_clusters.txt', 'w') as f:
    for line in remove_list:
        f.write(f"{line}\n")


rule all: 
    input: 
    	["exclude_datasets_no_clusters.txt"]
