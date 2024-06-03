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
parser = ArgumentParser(prog = 'inner-run-vespucci.R')
grid = read.delim("sh/grids/published_data/run-vespucci.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))
args = parser$parse_args()

source('R/functions/get_comparisons.R')

sc = readRDS(args$seurat_filename)
input = GetAssayData(sc, slot='counts')
meta = sc@meta.data

comparisons = get_comparisons(args$dataset, input, meta)

output_list = list()

go_list = ifelse(args$run_go, readRDS('data/metadata/go_list.rds'), NULL)
coords_col = ifelse(args$3d, c('x', 'y'), c('x', 'y', 'z'))
comparison_idx = 1
expr0 = comparisons[[comparison_idx]]$expr
meta0 = comparisons[[comparison_idx]]$meta

comparison = names(comparisons)[comparison_idx]
message("Processing comparison: ", comparison)

ves_res = run_vespucci(
	input = expr0,
	meta = meta0,
	go_list = go_list,
	coords_col = coords_col
)	

saveRDS(ves_res, args$output_filename)