setwd("~/git/vespucci-analysis")
options(stringsAsFactors = FALSE)
library(argparse)
library(Seurat)
library(tidyverse)
library(magrittr)
library(RANN)
library(Matrix)
library(sparseMatrixStats)
library(Vespucci)
# dynamically parse arguments

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-write-GO-matrix.R')
grid = read.delim("sh/grids/published_data/write-GO-matrix.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))

args = parser$parse_args()
print(args)

source('R/functions/get_comparisons.R')

sc = readRDS(args$seurat_filename)
expr = GetAssayData(sc, slot='counts')
meta = sc@meta.data

comparisons = get_comparisons(args$dataset, input, meta)

go_list = readRDS('data/metadata/go_list.rds')
go_mat = get_GO_matrix(comparisons[[1]]$expr, go_list)

new_sc = CreateSeuratObject(go_mat, meta=comparisons[[1]]$meta)
saveRDS(new_sc, args$output_filename)