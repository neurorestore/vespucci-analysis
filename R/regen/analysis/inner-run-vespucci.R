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
grid = read.delim("sh/grids/regen/run-vespucci.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))
args = parser$parse_args()

source('R/functions/get_comparisons.R')

sc = readRDS(args$seurat_filename)
input = GetAssayData(sc, slot='counts')
meta = sc@meta.data

output_list = list()

go_list = ifelse(args$run_go, readRDS('data/metadata/go_list.rds'), NULL)

ves_res = run_vespucci(
	input = input,
	meta = meta,
	go_list = go_list
)	

saveRDS(ves_res, args$output_filename)