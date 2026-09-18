Sys.setenv(
  RETICULATE_PYTHON = "/home/brouwer2/.conda/envs/R_CEPO/bin/python"
)
# script to create a gene ranking based on CEPO - writes values that still need to be ranked for gene-ranking purposes, and top 10% filtered to use in combination with LDSC (Li et al., medrxiv)
library(Seurat)
library(anndata)
library(Cepo)
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser = ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--infile", help ="Path to the input file .h5ad format.")
parser$add_argument("--dataset", help="DatasetID")
parser$add_argument("--level", help="level")
parser$add_argument("--outdir")

args = parser$parse_args()
h5ad=args$infile
level=args$level
outdir=args$outdir
dataset=args$dataset


ad = anndata::read_h5ad(h5ad)
# filter on cells (MT genes and/or too little transcripts - column generated in preprocessing
level_keep_column = sprintf("keep_%s",level)
#ad = ad[ad$obs[,level_keep_column]==TRUE,]


# cepo does not work if cell type clusters are smaller than 20, remove all cells of small/empty clusters
cell_type_column=sprintf("cell_type_%s",level)
ct=ad$obs[,cell_type_column]

tabs = table(ct)
print(tabs)
remove=names(tabs)[tabs<20]
ad = ad[!(ad$obs[,cell_type_column] %in% remove)] 
ct_filtered=ad$obs[,cell_type_column]
# remove genes that do not have ensembl ID
ad=ad[,!is.na(ad$var[,'ensembl'])]
rownames(ad$var)=ad$var[,'ensembl']


# log normalise CPM - this is a Seurat function that expects rows to be genes and columns to be cells - we need to transpose the h5ad matrix
out=NormalizeData(t(ad$X),scale.factor=1e6)

ds_res=Cepo(exprsMat=out,ct_filtered,prefilter_pzero = 0.4, exprsPct=0.05)
metrics=cbind(rownames(ds_res$stats),ds_res$stats)
colnames(metrics)[1] = "GENE"
colnames(metrics) <- gsub("\\.", "_", names(metrics))

write.table(metrics, sprintf("%s/%s.%s_cepo.tsv", outdir,dataset,level), quote = F, sep = "\t", row.names = F)

# Cepo+LDSC only includes genes that have > 1 TPM within the cell type - we create a filtered matrix with NAs for lowly expressed genes
# Cepo can only be computed for all cell types together, per gene, so this is done first - we cannot filter on genes a priori as this would effect the values in other cell types

out_scaled=NormalizeData(t(ad$X),scale.factor=1e6,normalization.method="RC")
countpM=aggregate(t(out_scaled),by=list(ct_filtered),FUN="mean")
rownames(countpM) = countpM[,1]
countpM=t(countpM[,2:dim(countpM)[2]]>0)
countpM_filtered = countpM[metrics[,1],]
countpM_filtered[countpM_filtered==FALSE] <- NA
countpM_filtered[countpM_filtered==TRUE] <- 1

metrics_filtered=as.matrix(metrics[,2:dim(metrics)[2]])*countpM_filtered
metrics_filtered = data.frame(cbind(metrics[,1],metrics_filtered))

colnames(metrics_filtered)[1] = "GENE"
write.table(metrics_filtered, sprintf("%s/%s.%s_cepo_filtered.tsv", outdir,dataset,level), quote = F, sep = "\t", row.names = F)

