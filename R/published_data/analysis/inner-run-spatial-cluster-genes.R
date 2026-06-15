# Run Magellan on each barcode within a given real dataset
setwd("~/git/vespucci-analysis/")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-run-spatial-cluster-genes-real-data.R')
grid = read.delim("sh/grids/rejected_review/run-spatial-cluster-genes-real-data.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))

args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(Seurat)
library(Vespucci)

sc = readRDS(args$seurat_file)
mat = GetAssayData(sc, slot='counts')
meta = sc@meta.data
ves_res = readRDS(args$ves_file)[[1]]

cluster_res = find_de_clusters(mat, input, ves_res)
saveRDS(cluster_res, args$output_filename)