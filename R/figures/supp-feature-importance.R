library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(ggpubr)

# setwd('C:/Users/teo/Documents/EPFL/projects/vespucci/')
setwd('~/git/vespucci-analysis/')
source('R/theme.R')

# comparison = 'young_old'
comparisons = c('young_old', 'treated_old')
for (comparison in comparisons) {
    # sc = readRDS(paste0('data/real_data/seurat/regen_final_', comparison, '.rds'))
    sc = readRDS(paste0('data/real_data/seurat/regen_final_', comparison, '.rds'))
    meta = sc@meta.data %>% mutate(barcode = gsub('-','_',barcode))
    ves_res = readRDS(paste0('data/real_data/vespucci/regen_final_', comparison, '-seed=42-with_feature_importance.rds'))[[1]]
    dat = ves_res$spatial_auc_result$feature_importance
    barcodes_with_auc = dat %>% pull(barcode) %>% unique()
    
    rostral_barcodes = meta %>% filter(x < -200) %>% pull(barcode)
    rostral_importance = dat %>% 
        filter(barcode %in% rostral_barcodes) %>%
        group_by(gene, barcode) %>%
        summarise(mean_importance = sum(mean_importance)/10) %>% # since there are 10 subsamples
        ungroup() %>%
        group_by(gene) %>%
        summarise(rostral_importance = mean(mean_importance)) %>% 
        ungroup()
    
    caudal_barcodes = meta %>% filter(x > 200) %>% pull(barcode)
    caudal_importance = dat %>% 
        filter(barcode %in% caudal_barcodes) %>%
        group_by(gene, barcode) %>%
        summarise(mean_importance = sum(mean_importance)/10) %>% # since there are 10 subsamples
        ungroup() %>%
        group_by(gene) %>%
        summarise(caudal_importance = mean(mean_importance)) %>% 
        ungroup()
    
    lesion_barcodes = meta %>% filter(x > -100 & x < 100) %>% pull(barcode)
    lesion_importance = dat %>% 
        filter(barcode %in% lesion_barcodes) %>%
        group_by(gene, barcode) %>%
        summarise(mean_importance = sum(mean_importance)/10) %>% # since there are 10 subsamples
        ungroup() %>%
        group_by(gene) %>%
        summarise(lesion_importance = mean(mean_importance)) %>% 
        ungroup()
    
    cor_df = rostral_importance %>% 
        full_join(caudal_importance, by='gene')
    cor_df[is.na(cor_df)] = 0
    
    p2 = cor_df %>%
        ggplot(aes(x=rostral_importance, y=caudal_importance)) +
        ggrastr::rasterise(
            geom_point(size = 0.1, shape='.'), dpi = 600
        ) +
        boxed_theme() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
#            axis.text.y = element_blank(),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(angle=45, hjust=1),
#            axis.ticks.x = element_blank(),
#            axis.ticks.y = element_blank(),
        ) +
        ggtitle('Rostral vs Caudal') +
        stat_cor(size=1.5, aes(label = ..r.label..))
    p2
    
    cor_df = rostral_importance %>% 
        full_join(lesion_importance, by='gene')
    cor_df[is.na(cor_df)] = 0
    
    p3 = cor_df %>%
        ggplot(aes(x=rostral_importance, y=lesion_importance)) +
        ggrastr::rasterise(
            geom_point(size = 0.1, shape='.'), dpi = 600
        ) +
        boxed_theme() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
#            axis.text.y = element_blank(),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(angle=45, hjust=1),
#            axis.ticks.x = element_blank(),
#            axis.ticks.y = element_blank(),
        ) +
        ggtitle('Rostral vs Lesion') +
        stat_cor(size=1.5, aes(label = ..r.label..))
    p3
    
    cor_df = caudal_importance %>% 
        full_join(lesion_importance, by='gene')
    cor_df[is.na(cor_df)] = 0
    
    p4 = cor_df %>%
        ggplot(aes(x=caudal_importance, y=lesion_importance)) +
        ggrastr::rasterise(
            geom_point(size = 0.1, shape='.'), dpi = 600
        ) +
        boxed_theme() +
        ggtitle('Caudal vs Lesion') +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
#            axis.text.y = element_blank(),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(angle=45, hjust=1),
#            axis.ticks.x = element_blank(),
#            axis.ticks.y = element_blank(),
        ) +
        stat_cor(size=1.5, aes(label = ..r.label..))
    p4
    
    # now against default augur feature importance
    augur_res = readRDS(paste0('data/real_data/regen/augur/regen_final_', comparison, '-augur_res.rds'))
    augur_importance = augur_res$feature_importance %>%
        group_by(gene) %>%
        summarise(augur_importance = mean(importance)) %>% 
        ungroup()
    
    cor_df = rostral_importance %>% 
        full_join(augur_importance, by='gene')
    cor_df[is.na(cor_df)] = 0
    p5 = cor_df %>%
        ggplot(aes(x=rostral_importance, y=augur_importance)) +
        ggrastr::rasterise(
            geom_point(size = 0.1, shape='.'), dpi = 600
        ) +
        boxed_theme() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
#            axis.text.y = element_blank(),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(angle=45, hjust=1),
#            axis.ticks.x = element_blank(),
#            axis.ticks.y = element_blank(),
        ) +
        ggtitle('Rostral vs Augur') +
        stat_cor(size=1.5, aes(label = ..r.label..))
    p5
    
    cor_df = caudal_importance %>% 
        full_join(augur_importance, by='gene')
    cor_df[is.na(cor_df)] = 0
    p6 = cor_df %>%
        ggplot(aes(x=caudal_importance, y=augur_importance)) +
        ggrastr::rasterise(
            geom_point(size = 0.1, shape='.'), dpi = 600
        ) +
        boxed_theme() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
#            axis.text.y = element_blank(),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(angle=45, hjust=1),
#            axis.ticks.x = element_blank(),
#            axis.ticks.y = element_blank(),
        ) +
        ggtitle('Caudual vs Augur') +
        stat_cor(size=1.5, aes(label = ..r.label..))
    p6
    
    cor_df = lesion_importance %>% 
        full_join(augur_importance, by='gene')
    cor_df[is.na(cor_df)] = 0
    p7 = cor_df %>%
        ggplot(aes(x=lesion_importance, y=augur_importance)) +
        ggrastr::rasterise(
            geom_point(size = 0.1, shape='.'), dpi = 600
        ) +
        boxed_theme() +
        ggtitle('Lesion vs Augur') +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
#            axis.text.y = element_blank(),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(angle=45, hjust=1),
#            axis.ticks.x = element_blank(),
#            axis.ticks.y = element_blank(),
        ) +
        stat_cor(size=1.5, aes(label = ..r.label..))
    p7
    
    # finally aggregate classifer
    aggregate_classifer_files = list.files('data/real_data/regen/aggregate_classifier/', full.names=T, pattern=comparison)
    agg_importance = map_df(aggregate_classifer_files, function(x){
        readRDS(x)$rf$importance %>% data.frame() %>% set_colnames('agg_importance') %>% rownames_to_column('gene')
    })
    agg_importance = agg_importance %>%
        group_by(gene) %>%
        summarise(agg_importance = mean(agg_importance)) %>% 
        ungroup()
    cor_df = rostral_importance %>% 
        full_join(agg_importance, by='gene')
    cor_df[is.na(cor_df)] = 0
    p8 = cor_df %>%
        ggplot(aes(x=rostral_importance, y=agg_importance)) +
        ggrastr::rasterise(
            geom_point(size = 0.1, shape='.'), dpi = 600
        ) +
        boxed_theme() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
#            axis.text.y = element_blank(),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(angle=45, hjust=1),
#            axis.ticks.x = element_blank(),
#            axis.ticks.y = element_blank(),
        ) +
        ggtitle('Rostral vs Random forest') +
        stat_cor(size=1.5, aes(label = ..r.label..))
    p8
    
    cor_df = caudal_importance %>% 
        full_join(agg_importance, by='gene')
    cor_df[is.na(cor_df)] = 0
    p9 = cor_df %>%
        ggplot(aes(x=caudal_importance, y=agg_importance)) +
        ggrastr::rasterise(
            geom_point(size = 0.1, shape='.'), dpi = 600
        ) +
        boxed_theme() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
#            axis.text.y = element_blank(),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(angle=45, hjust=1),
#            axis.ticks.x = element_blank(),
#            axis.ticks.y = element_blank(),
        ) +
        ggtitle('Caudual vs Random forest') +
        stat_cor(size=1.5, aes(label = ..r.label..))
    p9
    
    cor_df = lesion_importance %>% 
        full_join(agg_importance, by='gene')
    cor_df[is.na(cor_df)] = 0
    p10 = cor_df %>%
        ggplot(aes(x=lesion_importance, y=agg_importance)) +
        ggrastr::rasterise(geom_point(size = 0.1, shape='.'), dpi = 600) +
        boxed_theme() +
        ggtitle('Lesion vs Random forest') +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
#            axis.text.y = element_blank(),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(angle=45, hjust=1),
#            axis.ticks.x = element_blank(),
#            axis.ticks.y = element_blank(),
        ) +
        stat_cor(size=1.5, aes(label = ..r.label..))
    p10
    
    out_p = wrap_plots(
        p2, p3, p4,
        p5, p6, p7,
        p8, p9, p10,
        nrow = 3
    )
    # out_p
    ggsave(paste0('fig/EFig16/', comparison, '-regions_cor_plot.pdf'), out_p, width=9, height=9, units='cm')
}
