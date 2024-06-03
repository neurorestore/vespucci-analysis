setwd('~/git/vespucci-analysis/')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(countsplit)
source('R/theme.R')

new_v_sc = readRDS('data/regen/rnaseq/spatial/seurat/aligned_seurat.rds')
expr = new_v_sc@assays$RNA@layers$counts
colnames(expr) = new_v_sc@assays$RNA@cells@dimnames[[1]]
rownames(expr) = new_v_sc@assays$RNA@features@dimnames[[1]]
sc = CreateSeuratObject(expr, meta.data=new_v_sc@meta.data)
rm(new_v_sc)
gc()

sc@meta.data %<>% mutate(replicate = paste0(sample, '_', label))
sc@meta.data$barcode = colnames(sc)
sc[["n_umis"]] = colSums(expr)
sc[["n_genes"]] = colSums(expr > 0)
sc[["pct_mito"]] <- PercentageFeatureSet(sc, pattern = "^mt-")

sc %<>%
    subset(
        replicate %in% c(
            '1_B_young',
            '2_B_young',
            '3_B_young',
            '4_A_young',
            '1_C_old',
            '2_D_old',
            '3_A_old',
            '4_D_old',
            '1_B_treated',
            '4_C_treated'
        )
    )
sc@meta.data %<>%
    mutate(
        x = ifelse(replicate == '2_D_treated', x - 60, x)
    )

meta = sc@meta.data
meta %<>% filter(
    x > -400,
    x < 400,
    y > -150,
    y < 150
)
sc = sc[,meta$barcode]
meta1 = sc@meta.data %>%
    filter(replicate %in% final_reps_to_keep)
sc1 = sc[,meta1$barcode]
saveRDS(sc1, 'data/regen/seurat/regen_final.rds')
