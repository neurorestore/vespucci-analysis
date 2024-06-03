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

sc = readRDS(args$seurat_filename)
go_list = readRDS('data/metadata/go_list.rds')
expr = GetAssayData(sc, slot='counts'
meta = sc@meta.data

meta_tmp = meta %>%
	filter(label %in% c('young', 'old'))
go_mat = get_GO_matrix(expr[,rownames(meta_tmp)], go_list)
new_sc = CreateSeuratObject(go_mat, meta=meta_tmp)
saveRDS(new_sc, 'data/regen/seurat_GO/regen_final_young_old.rds')

meta_tmp = meta %>%
	filter(label %in% c('treated', 'old'))
go_mat = get_GO_matrix(expr[,rownames(meta_tmp)], go_list)
new_sc = CreateSeuratObject(go_mat, meta=meta_tmp)
saveRDS(new_sc, 'data/regen/seurat_GO/regen_final_treated_old.rds')