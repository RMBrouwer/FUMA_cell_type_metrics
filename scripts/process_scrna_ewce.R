Sys.setenv(
  RETICULATE_PYTHON = "/home/brouwer2/.conda/envs/R_with_python/bin/python"
)

library(EWCE)
library(anndata)
suppressPackageStartupMessages(library("argparse"))

set.seed(1234)

# create parser object
parser = ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--h5ad", help ="Path to the h5ad file.")
parser$add_argument("--outdir", help ="Path to the output file")
parser$add_argument("--file_prefix", help ="The prefix of the output file") #for example: 1_Allen_MCA_Human_2019
parser$add_argument("--level", help="Cell type level, e.g. 'level_1' (assuming the cell types are labeled 'cell_type_level_1' etc in the dataset")
args = parser$parse_args()

# get arguments
h5ad = args$h5ad
outdir = args$outdir
file_prefix = args$file_prefix
level = args$level


# get count from anndata
adata = read_h5ad(h5ad)
rownames(adata$var) = adata$var$ensembl

#remove all cells of small/empty clusters < 20
cell_type_column=sprintf("cell_type_%s",level)
ct=adata$obs[,cell_type_column]

tabs = table(ct)
print(tabs)
remove=names(tabs)[tabs<20]
adata = adata[!(adata$obs[,cell_type_column] %in% remove)]

# remove genes that do not have ensembl ID
adata=adata[,!is.na(adata$var[,'ensembl'])]
rownames(adata$var)=adata$var[,'ensembl']



# Turn counts into matrix
# Because as.matrix(adata$X) returns a matrix where rows are the cells and columns are the genes. Therefore, we need to transpose it to have rows are the cells and columns are the genes.
count <- t(as.matrix(adata$X))
 print(ncol(count))
 print(nrow(count))

if (level == "level_1")
	annotLevels = list(level1class=adata$obs$cell_type_level_1)
if (level == "level_2")
	annotLevels = list(level1class=adata$obs$cell_type_level_2)
if (level == "level_3")
        annotLevels = list(level1class=adata$obs$cell_type_level_3)
annotLevels = list(level1class=adata$obs[sprintf("cell_type_%s",level)][,1])

specificity_out <- EWCE::generate_celltype_data(
  exp = count,
  annotLevels = annotLevels,
  groupName = level,
  savePath = outdir,
  file_prefix = file_prefix,
  input_species = "human",
  dendrograms = FALSE
)
