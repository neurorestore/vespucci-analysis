# use project/envs/nebula
library(tidyverse)
library(magrittr)
library(Seurat)
library(hdf5r)

outer_dir = 'data/regen/rnaseq/aligned/'
samples = basename(list.dirs(outer_dir, recursive=F))

sc = NULL
for (sample in samples) {
	print(sample)	
	data_dir = paste0(outer_dir, sample, '/', sample, '/outs/')
	sc0 = Load10X_Spatial(data_dir)
	sc0 %<>% RenameCells(sample)
	if (is.null(sc)) {
		sc = sc0
	} else {
		sc = merge(sc, sc0)
	}
}

saveRDS(sc, 'data/regen/rnaseq/spatial/seurat/raw_seurat.rds')