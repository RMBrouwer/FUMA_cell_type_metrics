# This script aims to get the output from ewce into a text file


library(EWCE)
library(data.table)
library(dplyr)
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser = ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--infile", help ="Path to the output file of the process_scrna_ewce.R script - rda format.")
parser$add_argument("--dataset", help="DatasetID")
parser$add_argument("--level", help="level")
parser$add_argument("--outdir")
args = parser$parse_args()


ctd=args$infile
level=args$level
outdir=args$outdir
dataset=args$dataset

ctd = EWCE::load_rdata(ctd)

# specificity
specificity_df = as.data.frame(as.matrix(ctd[[1]]$specificity))

#replace space with _
names(specificity_df) <- gsub(" ", "_", names(specificity_df))

#convert rowname to a column
# specificity_df$symbol <- row.names(specificity_df)
specificity_df$GENE <- row.names(specificity_df)

# we do not need the following since names are already in ensembl ID - but this is not stored in the rda file. Will redo
# convert gene symbol to ensemble
#conversion_file = fread("/projects/0/vusr0480/Preprocessing_scRNA/code/conversion_files/gene_names_human.txt", quote="")

#conversion_file_short = conversion_file %>%
#  select("Approved symbol", "Ensembl gene ID")
#colnames(conversion_file_short) = c("symbol", "GENE")

#specificity_df_ensembl = merge(specificity_df, conversion_file_short, by = c("symbol"))
#colnames(specificity_df_ensembl)[which(colnames(specificity_df_ensembl) == "ensembl")] = "GENE"

specificity_df_ensembl_fmt = specificity_df %>%
#  select(-symbol) %>%
  select(GENE, everything())
#write.table(specificity_df, sprintf("%s/%s.%s_spec_ewce_linear.tsv", outdir,dataset,level), quote = F, sep = "\t", row.names = FALSE)
write.table(specificity_df_ensembl_fmt, sprintf("%s/%s.%s_spec_ewce_linear.tsv", outdir,dataset,level), quote = F, sep = "\t", row.names = FALSE)

# mean
mean_df = as.data.frame(as.matrix(ctd[[1]]$mean_exp))

#replace space with _
names(mean_df) <- gsub(" ", "_", names(mean_df))

#convert rowname to a column
#mean_df$symbol <- row.names(mean_df)
mean_df$GENE <- row.names(mean_df)

#mean_df_ensembl = merge(mean_df, conversion_file_short, by = c("symbol"))
#colnames(mean_df_ensembl)[which(colnames(mean_df_ensembl) == "ensembl")] = "GENE"

mean_df_ensembl_fmt = mean_df %>%
#  select(-symbol) %>%
  select(GENE, everything())

mean_df_ensembl_fmt$Average = rowMeans(mean_df_ensembl_fmt[,2:dim(mean_df_ensembl_fmt)[2]])

write.table(mean_df_ensembl_fmt, sprintf("%s/%s.%s_mean_ewce_linear.tsv", outdir,dataset,level), quote = F, sep = "\t", row.names = FALSE)


