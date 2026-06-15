setwd('~/git/vespucci-analysis/')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(sparseMatrixStats)

source('R/theme.R')
source('R/functions/utils.R')

# just old treated
go_res = readRDS('data/real_data/spatial_cluster_genes/gene_set_enrichment/regen_final_treated_old.rds')
meta0 = readRDS('data/real_data/meta/regen_final.rds')

label = 'treated'
cluster_res = readRDS(paste0('data/real_data/spatial_cluster_genes/regen_final_treated_old-norm=norm-clust=kmeans-label=', label, '-seed=42-expr=module.rds'))
barcode_expr_df = cluster_res$barcode_expr_df
genes_cluster_df = cluster_res$genes_cluster_df
plot_list = list()
for (cluster in sort(unique(barcode_expr_df$cluster))) {
    genes = genes_cluster_df %>% filter(cluster == !!cluster) %>% pull(gene)
    gos = go_res %>% filter(gsub('cluster_', '', cluster) == !!cluster, label == !!label) %>% arrange(p_val_adj) %>% head(10) %>% arrange(p_val)
    gos %<>% mutate(go_name = factor(go_name, levels=rev(gos$go_name)))
    meta = meta0 %>% inner_join(barcode_expr_df %>% filter(cluster == !!cluster))        
    ## by condition
    conditions = unique(meta$label)
    coords = map_dfr(conditions, ~ {
        condition = .x
        fit = loess(expr ~ x * y, data = filter(meta, label == condition), span = 0.02, degree = 1)
        coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5)) %>%
            mutate(interp = predict(fit, newdata = .)) %>%
            drop_na(interp) %>%
            mutate(condition = condition)
    })

    coords %<>% mutate(condition = factor(str_to_title(condition), levels=c('Young', 'Treated', 'Old')))
    # range = range(coords$interp)
    range = quantile(coords$interp, c(0.01, 0.99))
    brks = c(range[1] + 0.1 * diff(range), range[2] - 0.1 * diff(range))
    labels = c('min', 'max')
    p1a = coords %>%
        mutate(interp = winsorize(interp, quantile(coords$interp, c(0.01, 0.99)))) %>%
        ggplot(aes(x = x, y = y, fill = interp)) +
        facet_grid(~ condition, switch='y') +
        geom_raster() +
        scale_y_reverse(expand = c(0, 0), breaks = c(-150, 150), labels = c(expression('R'%->%""),expression(""%<-%'L'))
        ) +
        scale_x_continuous(expand = c(0, 0), breaks = c(-400, 400),labels = c(expression(""%<-%'Rostral'),expression('Caudal'%->%""))
        ) +
        scale_fill_gradientn(name = 'Expression', colors = nr_tree_red, limits = range, breaks = brks, labels = labels) +
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
        coord_fixed() +
        boxed_theme() +
        ggtitle(paste0('Cluster ', cluster, '\n# of genes: ', length(genes))) +
        theme(
            plot.title = element_text(size=4),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(angle = 90, hjust = c(1, 0)),
            axis.text.x = element_text(hjust = c(0, 1)),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.length.x = unit(0, 'lines'),
            axis.ticks.length.y = unit(0, 'lines'),
            legend.position = 'none',
            legend.justification = 'bottom',
            legend.key.width = unit(0.12, 'lines'),
            legend.key.height = unit(0.12, 'lines')
        )
    p1b = gos %>% 
        ggplot(aes(x = go_name, y = -log10(p_val))) +
        geom_hline(aes(yintercept = 0), size = 0.4, color = 'grey88') +
        geom_segment(aes(xend = go_name, yend = 0), color = 'grey88') +
        scale_x_discrete(position='top')+
        geom_point(shape = 21, stroke = 0.2, size = 0.9, color = 'black', fill = 'grey80') +
        geom_segment(aes(y = -Inf, yend = Inf, x = 10.5, xend = 10.5),  color = 'grey20', size = 0.225) +
        scale_y_continuous(expression(-log[10]~'p'), expand=c(0, 5)) +
        guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
        coord_flip() +
        boxed_theme() +
        theme(
            # strip.text = element_blank(),
            aspect.ratio = 1.8,
            axis.title.y = element_blank(),
            legend.position = 'right',
            legend.justification = 'bottom',
            legend.key.width = unit(0.2, 'lines'),
            legend.key.height = unit(0.2, 'lines'),
            legend.text = element_text(size = 5),
            legend.title = element_text(size = 5),
        )
    plot_list[[length(plot_list)+1]] = wrap_plots(p1a,p1b,nrow=1,widths=c(1,1.5))
}
out_p = wrap_plots(plot_list, ncol=1)
ggsave(paste0('fig/EFig14/regen_final_treated_old-label=', label, '.pdf'), out_p, width = 15, height = 30, units = 'cm', useDingbats = FALSE)

# first correlation between each cluster
run_grid = go_res %>% dplyr::select(label, cluster) %>% distinct() %>% dplyr::rename(cluster1 = cluster) %>% left_join(go_res %>% dplyr::select(label, cluster) %>% distinct() %>% dplyr::rename(cluster2 = cluster), relationship = 'many-to-many') %>% filter(cluster1 != cluster2)

cor_df = map_df(1:nrow(run_grid), function(i){
    tmp_row = run_grid[i,]
    tmp_df = go_res %>% filter(label == tmp_row$label, cluster == tmp_row$cluster1) %>% dplyr::select(go_name, p_val) %>% dplyr::rename(cluster1_pval = p_val) %>% left_join(go_res %>% filter(label == tmp_row$label, cluster == tmp_row$cluster2) %>% dplyr::select(go_name, p_val) %>% dplyr::rename(cluster2_pval = p_val))
    tmp_row %>% 
        mutate(
            pearson_cor = cor(-log(tmp_df$cluster1_pval, 10), -log(tmp_df$cluster2_pval, 10), method='pearson'),
            spearman_cor = cor(-log(tmp_df$cluster1_pval, 10), -log(tmp_df$cluster2_pval, 10), method='spearman')
        )
})

plot_list = list()
for (label in c('old', 'treated')) {
    cor_df1 = cor_df %>% filter(label == !!label) %>% mutate(cluster1 = gsub('cluster_', 'Cluster ', cluster1), cluster2 = gsub('cluster_', 'Cluster ', cluster2))
    # range = range(cor_df1$pearson_cor)
    range = c(-1,1)
    # range = range(cor_df1$spearman_cor)
    brks = range

    p1 = cor_df1 %>%
        ggplot(aes(x = cluster1, y = cluster2)) +
        geom_tile(color = 'white', aes(fill = pearson_cor)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",name = 'Corr.', limits=range,breaks = brks, labels = c(-1,1)) +
        guides(fill = guide_colorbar(theme = theme(legend.ticks = element_blank(),legend.frame = element_rect(colour = "black", size = 0.25))), color = guide_colorbar(theme = theme(legend.ticks = element_blank(),legend.frame = element_rect(colour = "black", size = 0.25)))) +
        coord_fixed() +
        boxed_theme() +
        ggtitle(str_to_title(label)) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.key.width = unit(0.18, 'lines'),
              legend.key.height = unit(0.15, 'lines'),
              legend.position = 'bottom',
              legend.justification = 'right')
    plot_list[[length(plot_list)+1]] = p1
}
out_p = wrap_plots(plot_list, nrow=1)
ggsave('fig/EFig14/gene_set_enrichment_pval_cor.pdf', out_p, width=6.5, height=6.5, units='cm')
