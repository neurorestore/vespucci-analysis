library(tidyverse)
library(magrittr)
library(Matrix)
library(Seurat)

source('R/functions/registration_functions.R')

sc = readRDS('published_data/raw_data/Kathe2022/Kathe2022_raw_seurat.rds')
meta = sc@meta.data
rownames(meta) = gsub('-','_',rownames(meta))
meta$barcode = rownames(meta)

meta$ori_x = sc@images$slice1@coordinates$imagecol
meta$ori_y = sc@images$slice1@coordinates$imagerow

meta %<>%
	mutate(
		label = group,
		replicate = mouse_id,
		x = abs(ori_x - 250),
		y = -ori_y
	) %>%
	filter(lowres == 'grey') %>%
	dplyr::select(ori_x, ori_y, x,y,label,replicate, barcode,lowres)

meta %<>%
	mutate(label = recode(label, 
		'6wNT' = 'SCI',
		'6wT' = 'EES_REHAB'
		)) %>%
	filter(label %in% c('EES_REHAB', 'SCI'))

counts = GetAssayData(sc, slot='counts')
colnames(counts) = gsub('-', '_', colnames(counts))
counts = counts[,rownames(meta)]
new_sc = CreateSeuratObject(counts, meta.data = meta)

saveRDS(new_sc, 'data/published_data/seurat/Kathe2022.rds')