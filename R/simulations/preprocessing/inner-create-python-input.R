setwd("~/git/vespucci-analysis")
options(stringsAsFactors = FALSE)
library(argparse)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(magrittr)
library(reticulate)
library(Matrix)

# dynamically parse arguments

parser = ArgumentParser(prog = 'inner-create-python-input.R')
grid = read.delim("sh/grids/simulations/create-python-input.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))

args = parser$parse_args()
print(args)

sc = readRDS(args$input_filename)
sc@assays$magellan = NULL
sc@assays$validation = NULL

# for sepal to work
sc@assays$originalexp@meta.features$name =  sc@assays$originalexp@meta.features$Gene

# save as AnnData 
anndata_dir = paste0(args$output_dir, 'anndata/')
if (!dir.exists(anndata_dir)) dir.create(anndata_dir, recursive=T)

h5_seurat_filename = paste0(anndata_dir, args$output_prefix, '.h5Seurat')
SaveH5Seurat(sc, filename = h5_seurat_filename)
Convert(h5_seurat_filename, dest = "h5ad")

# save as mtx
mtx_dir = args$mtx_output_dir
if (!dir.exists(mtx_dir)) dir.create(mtx_dir, recursive=T)
counts = GetAssayData(sc, slot='counts')
writeMM(counts, args$mtx_output_filename)

meta_dir = paste0(args$output_dir, 'meta/')
if (!dir.exists(meta_dir)) dir.create(meta_dir, recursive=T)
meta = sc@meta.data
meta_output_filename = paste0(meta_dir, args$output_prefix, '_metadata.csv.gz')
write.csv(meta, gzfile(meta_output_filename))
coords_only = meta %>% select(x,y)
coords_output_filename = paste0(meta_dir, args$output_prefix, '_coords.csv.gz')
write.csv(coords_only, gzfile(coords_output_filename))