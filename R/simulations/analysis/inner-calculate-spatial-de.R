setwd("~/git/vespucci-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-calculate-spatial-de.R')
grid = read.delim("sh/grids/simulations/calculate-spatial-de.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))

args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
source('R/functions/spatial_de_methods.R')

sc = readRDS(args$input_filename)

input = GetAssayData(sc, slot='counts')
input = input[rowSums(input) > 0,]
meta = sc@meta.data
gene_features = sc@assays$originalexp@meta.features

selected_features = NULL
split_size = as.integer(str_split(args$split_size, '_')[[1]][2])
if (args$split_size != '1') {	
	split_size_index = as.integer(str_split(args$split_size, '_')[[1]][1])
	features_chunks = split(rownames(input), cut(seq_along(rownames(input)), split_size, labels = FALSE))
	selected_features = features_chunks[[split_size_index]]
}

giotto_python_path = 'envs/de-methods/bin/python'
output_res = run_spatial_de(input, meta, gene_features, de_method = args$de_method, selected_features = selected_features, giotto_python_path = giotto_python_path)
saveRDS(output_res, args$output_filename)