setwd("~/git/vespucci")
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
parser = ArgumentParser(prog = 'inner-run-vespucci-feature-importance.R')
grid = read.delim("sh/grids/rejected_review/run-vespucci-feature-importance.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))

args = parser$parse_args()
print(args)

source('R/functions/get_comparisons.R')
sc = readRDS(args$seurat_file)
run_countsplit = T
input = GetAssayData(sc, slot='counts')
input[is.na(input)] = 0

meta = sc@meta.data
# filter zero counts barcodes (due to count splitting)
input = input[,colSums(input) > 0]
input = input[,colnames(input) %in% rownames(meta)]
meta = meta[colnames(input),]

comparisons = get_comparisons(args$dataset, input, meta)
rm(sc, input, meta)
gc()

output_list = list()

for (comparison_idx in seq_along(comparisons)) {
	expr0 = comparisons[[comparison_idx]]$expr
	colnames(expr0) = gsub('-', '_', colnames(expr0))

	meta0 = comparisons[[comparison_idx]]$meta
	rownames(meta0) = gsub('-', '_', rownames(meta0))
	meta0$barcode = rownames(meta0)

	comparison = names(comparisons)[comparison_idx]
	message("Processing comparison: ", comparison)

	expr0 = expr0[,meta0$barcode]

	if (grepl('3d', args$dataset)) {
		coord_cols = c('x','y','z')
	} else {
		coord_cols = c('x','y')
	}
	options(future.globals.maxSize = 8000 * 1024^2)
	ves_res = run_vespucci(
		expr0, 
		meta0,
		coord_cols = coord_cols,
		seed = args$seed,
		save_feature_importance = T
	)

	output_list[[comparison]] = ves_res
}	

saveRDS(output_list, args$output_filename)