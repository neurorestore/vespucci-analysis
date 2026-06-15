library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(cetcolor)
library(pROC)
library(PRROC)

source('R/theme.R')

# test different seeds
dat = readRDS('data/real_data/summaries/different_seeds/real_data-auc-summary.rds') %>%
    filter(dataset %in% c('Calcagno2022', 'regen_final_treated_old', 'regen_final_young_old')) %>% 
    mutate(
        cor = pearson_cor,
        dataset = case_when(
            dataset == 'regen_final_treated_old' ~ 'Treated vs old SCI',
            dataset == 'regen_final_young_old' ~ 'Young vs old SCI',
            T ~ dataset
        )
    )

labs = dat %>%
    group_by(dataset) %>%
    summarize(
        stats_val = mean(cor),
        val = max(cor)
    ) %>% 
    ungroup() %>%
    mutate(
        text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 2)),
        text_y = ifelse(val < 0, 0, val)
    )
pal = nr_base_3 %>% setNames(unique(dat$dataset))

p1 = dat %>%
    ggplot(aes(x = dataset, y = cor)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA, fill='grey60') +
    geom_label(dat = labs,
               aes(label = text_val, y = text_y), 
               color = ifelse(labs$text_val == 'OOT', 'grey', 'black'),
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0, vjust = 0.5,
               show.legend = FALSE) +
    coord_flip() +
    ggtitle(expression('AUCs correlation')) +
    scale_y_continuous(breaks = seq(0.8, 1, 0.1), limits = c(0.8, 1.05)) +
    boxed_theme() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        aspect.ratio = 0.45,
        legend.position = 'none',
        legend.justification = 'bottom',
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.15, 'lines')
    )
p1


dat = readRDS('data/real_data/summaries/different_seeds/real_data-pvals-summary.rds') %>%
    filter(dataset %in% c('Calcagno2022', 'regen_final_treated_old', 'regen_final_young_old')) %>% 
    mutate(
        cor = spearman_cor,
        dataset = case_when(
            dataset == 'regen_final_treated_old' ~ 'Treated vs old SCI',
            dataset == 'regen_final_young_old' ~ 'Young vs old SCI',
            T ~ dataset
        )
    )
labs = dat %>%
    group_by(dataset) %>%
    summarize(
        stats_val = mean(cor),
        val = max(cor)
    ) %>% 
    ungroup() %>%
    mutate(
        text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 2)),
        text_y = ifelse(val < 0, 0, val)
    )
pal = nr_base_3 %>% setNames(unique(dat$dataset))

p2 = dat %>%
    ggplot(aes(x = dataset, y = cor)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA, fill='grey60') +
    geom_label(dat = labs,
               aes(label = text_val, y = text_y), 
               color = ifelse(labs$text_val == 'OOT', 'grey', 'black'),
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0, vjust = 0.5,
               show.legend = FALSE) +
    coord_flip() +
    ggtitle(expression('P-values correlation')) +
    scale_y_continuous(breaks = seq(0.7, 1, 0.1), limits = c(0.7, 1.05)) +
    boxed_theme() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        aspect.ratio = 0.45,
        legend.position = 'none',
        legend.justification = 'bottom',
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.15, 'lines')
    )
p2
out_p = wrap_plots(p1, p2, nrow=2)
ggsave('fig/EFig6/real-data-different_seeds.pdf', out_p, width=4, height=6, units='cm')


datasets = c('Calcagno2022', 'regen_final_treated_old', 'regen_final_young_old')
max_cells = c(10, 50, 100, 200, 500, 1000)
for (dataset in datasets) {
    if (!grepl('regen_final', dataset)) {
        meta = readRDS(paste0('data/real_data/meta/', dataset, '.rds'))    
    } else {
        meta = readRDS('data/real_data/meta/regen_final.rds')
    }
    plot_list = list()
    for (curr_max_cells in max_cells) {
        if (curr_max_cells != 1000) {
            ves_res = readRDS(paste0('data/real_data/vespucci/', dataset, '-seed=42-nsub=10-bc=', curr_max_cells, '.rds'))[[1]]$spatial_auc_result$aucs
        } else {
            ves_res = readRDS(paste0('data/real_data/vespucci/', dataset, '-seed=42-nsub=10.rds'))[[1]]$spatial_auc_result$aucs
        }
        dat0 = meta %>% dplyr::select(barcode, label, replicate, x, y) %>% mutate(barcode = gsub('-', '_', barcode)) %>% inner_join(ves_res)
        if (grepl('regen_final', dataset)) {
            fixed_coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5))
            fit = loess(auc ~ x * y, data = dat0, span = 0.02, degree = 1)
            
            dat1 = fixed_coords %>%
                mutate(
                    auc_fit = predict(fit, .)
                ) %>%
                filter(!is.na(auc_fit))
            
            range = range(dat1$auc_fit)
            brks = c(range[1] + 0.1 * diff(range),
                     range[2] - 0.1 * diff(range))
            labels = format(range, digits = 2)
            new_color_pal = cet_pal(100, name = 'l19') %>% rev()
            
            p3_0 = dat1 %>%
                arrange(auc_fit) %>%
                ggplot(aes(x = x, y = y, fill = auc_fit, color = auc_fit)) +
                ggrastr::rasterise(
                    geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
                ) +
                scale_y_continuous(expand = c(0,0)) +
                scale_x_continuous(expand = c(0,0)) +
                scale_fill_gradientn(colours = new_color_pal,
                                     name = 'AUC', labels = labels,
                                     limits = range, breaks = brks) +
                scale_color_gradientn(colours = new_color_pal,
                                      name = 'AUC', labels = labels,
                                      limits = range, breaks = brks) +
                guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE, title.position = 'top'),
                       color = guide_colorbar(frame.colour = 'black', ticks = FALSE, title.position = 'top')) +
                coord_fixed() +
                boxed_theme() +
                ggtitle(paste0('Max barcodes: ', curr_max_cells)) +
                # geom_vline(xintercept = 0, size=0.1, linetype='dashed') +
                # geom_hline(yintercept = 0, size=0.1, linetype='dashed') +
                theme(
                    # aspect.ratio = 0.5,
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    # axis.text.y = element_text(angle = 90, hjust = c(1, 0)),
                    # axis.text.x = element_text(hjust = c(0, 1)),
                    axis.text.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.ticks.length.x = unit(0, 'lines'),
                    axis.ticks.length.y = unit(0, 'lines'),
                    legend.position = 'bottom',
                    legend.justification = 'right',
                    legend.key.width = unit(0.18, 'lines'),
                    legend.key.height = unit(0.18, 'lines'),
                    plot.title = element_text(size = 5)
                )
            p3_0
            plot_list[[length(plot_list)+1]] = p3_0
        } else {
            # interpolate in 2D
            fit = loess(auc ~ x * y, data = dat0, span = 0.015)
            dat0$auc_fit = predict(fit, dat0)
            
            range = range(dat0$auc_fit)
            brks = c(range[1] + 0.1 * diff(range),
                     range[2] - 0.1 * diff(range))
            labels = format(range, digits = 2)
            labels = c(paste0(labels[1], ' '),
                       paste0(' ', labels[2]))
            
            auc_pal = cet_pal(100, name = 'l19') %>% rev()
            
            p3_0 = dat0 %>%
                # arrange(auc_fit) %>%
                ggplot(aes(x = x, y = y, fill = auc_fit, color = auc_fit)) +
                ggtitle('Vespucci') +
                ggrastr::rasterise(
                    geom_point(size = 0.3, shape = 21, stroke = 0), dpi = 600
                ) +
                scale_y_continuous(expand = c(0, 0)) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_fill_gradientn(colours = auc_pal,
                                     name = 'AUC   ', labels = labels,
                                     limits = range, breaks = brks) +
                scale_color_gradientn(colours = auc_pal,
                                      name = 'AUC   ', labels = labels,
                                      limits = range, breaks = brks) +
                guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE, title.position = 'top'),
                       color = guide_colorbar(frame.colour = 'black', ticks = FALSE, title.position = 'top')) +
                ggtitle(paste0('Max barcodes: ', curr_max_cells)) +
                # coord_fixed() +
                boxed_theme(size_lg = 5, size_sm = 5) +
                theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      # axis.text.y = element_text(angle = 90, hjust = c(1, 0)),
                      # axis.text.x = element_text(hjust = c(0, 1)),
                      axis.text.y = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.ticks.length.x = unit(0, 'lines'),
                      axis.ticks.length.y = unit(0, 'lines'),
                      legend.position = 'bottom',
                      legend.justification = 'right',
                      legend.key.width = unit(0.18, 'lines'),
                      legend.key.height = unit(0.18, 'lines')
                )
            plot_list[[length(plot_list)+1]] = p3_0
        }
    }
    out_p = wrap_plots(plot_list, nrow=2)
    if (grepl('regen_final', dataset)) {
        ggsave(paste0('fig/EFig6/different_subsample_setup_', dataset, '.pdf'), out_p, width=7, height=6.5, units='cm')    
    } else {
        ggsave(paste0('fig/EFig6/different_subsample_setup_', dataset, '.pdf'), out_p, width=6, height=6.5, units='cm')
    }
}

# plot auc correlation grid
run_grid = tidyr::crossing(
    dataset = c('Calcagno2022', 'regen_final_treated_old', 'regen_final_young_old'),
    max_cells1 = c(10, 50, 100, 200, 500, 1000),
    max_cells2 = c(10, 50, 100, 200, 500, 1000)
) %>% filter(max_cells1 < max_cells2)

dat2 = map_df(1:nrow(run_grid), function(i){
    tmp_row = run_grid[i,]
    ves_res1 = readRDS(paste0('data/real_data/vespucci/', tmp_row$dataset, '-seed=42-nsub=10-bc=', tmp_row$max_cells1, '.rds'))[[1]]
    if (tmp_row$max_cells2 == 1000) {
        ves_res2 = readRDS(paste0('data/real_data/vespucci/', tmp_row$dataset, '-seed=42-nsub=10.rds'))[[1]]
    } else {
        ves_res2 = readRDS(paste0('data/real_data/vespucci/', tmp_row$dataset, '-seed=42-nsub=10-bc=', tmp_row$max_cells2, '.rds'))[[1]]
    }
    combined_auc_df = ves_res1$spatial_auc_result$aucs %>% dplyr::rename(auc1 = auc) %>% inner_join(ves_res2$spatial_auc_result$aucs %>% dplyr::rename(auc2 = auc), by='barcode')
    combined_pval_df = ves_res1$de_feature_result %>% dplyr::rename(p_val1 = p_val) %>% dplyr::select(feature, p_val1) %>% inner_join(ves_res2$de_feature_result %>% dplyr::rename(p_val2 = p_val) %>% dplyr::select(feature, p_val2), by='feature')
    tmp_row %>% mutate(
        auc_cor = cor(combined_auc_df$auc1, combined_auc_df$auc2, method='pearson', use='complete.obs'),
        pval_cor = cor(combined_pval_df$p_val1, combined_pval_df$p_val2, method='spearman', use='complete.obs')
    )
})

max_cells = c(10, 50, 100, 200, 500, 1000)
plot_list = list()
for (dataset in c('Calcagno2022', 'regen_final_treated_old', 'regen_final_young_old')) {
    title_name = case_when(
        dataset == 'regen_final_treated_old' ~ 'Treated vs. old',
        dataset == 'regen_final_young_old' ~ 'Young vs. old',
        T ~ dataset
    )
    dat2_1 = rbind(dat2 %>% filter(dataset == !!dataset), dat2 %>% filter(dataset == !!dataset) %>% set_colnames(c('dataset', 'max_cells2', 'max_cells1', 'auc_cor', 'pval_cor')) %>% dplyr::select(c('dataset', 'max_cells1', 'max_cells2', 'auc_cor', 'pval_cor'))) %>% mutate(val = auc_cor, max_cells1 = factor(max_cells1, levels=max_cells),max_cells2 = factor(max_cells2, levels=max_cells))
    range = range(dat2_1$val)
    brks = c(range[1] + 0.1 * diff(range),
             range[2] - 0.1 * diff(range))
    
    p4_0 = dat2_1 %>% 
        ggplot(aes(x = max_cells1, y = max_cells2)) +
        geom_tile(color = 'white', aes(fill = val)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",name = 'Corr.',breaks = brks, labels = format(range, digits = 2)) +
        guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black', title.position = 'top')) +
        coord_fixed() +
        boxed_theme() +
        ggtitle(title_name) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.key.width = unit(0.18, 'lines'),
              legend.key.height = unit(0.15, 'lines'),
              legend.position = 'bottom',
              legend.justification = 'right')
    plot_list[[length(plot_list)+1]] = p4_0
}
out_p = wrap_plots(plot_list, nrow=1)
ggsave('fig/EFig6/different_subsample_auc_cor.pdf', out_p, width=7, height=6.5, units='cm')

plot_list = list()
for (dataset in c('Calcagno2022', 'regen_final_treated_old', 'regen_final_young_old')) {
    title_name = case_when(
        dataset == 'regen_final_treated_old' ~ 'Treated vs. old',
        dataset == 'regen_final_young_old' ~ 'Young vs. old',
        T ~ dataset
    )
    dat2_1 = rbind(dat2 %>% filter(dataset == !!dataset), dat2 %>% filter(dataset == !!dataset) %>% set_colnames(c('dataset', 'max_cells2', 'max_cells1', 'auc_cor', 'pval_cor')) %>% dplyr::select(c('dataset', 'max_cells1', 'max_cells2', 'auc_cor', 'pval_cor'))) %>% mutate(val = pval_cor, max_cells1 = factor(max_cells1, levels=max_cells),max_cells2 = factor(max_cells2, levels=max_cells))
    range = range(dat2_1$val)
    brks = c(range[1] + 0.1 * diff(range),
             range[2] - 0.1 * diff(range))
    
    p4_0 = dat2_1 %>% 
        ggplot(aes(x = max_cells1, y = max_cells2)) +
        geom_tile(color = 'white', aes(fill = val)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",name = 'Corr.',breaks = brks, labels = format(range, digits = 2)) +
        guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black', title.position = 'top')) +
        coord_fixed() +
        boxed_theme() +
        ggtitle(title_name) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.key.width = unit(0.18, 'lines'),
              legend.key.height = unit(0.15, 'lines'),
              legend.position = 'bottom',
              legend.justification = 'right')
    p4_0
    plot_list[[length(plot_list)+1]] = p4_0
}
out_p = wrap_plots(plot_list, nrow=1)
ggsave('fig/EFig6/different_subsample_pval_cor.pdf', out_p, width=7, height=6.5, units='cm')
