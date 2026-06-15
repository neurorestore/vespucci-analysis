library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(cetcolor)
# library(pROC)
# library(PRROC)
library(ggpubr)

# setwd('C:/Users/teo/Documents/EPFL/projects/vespucci/')
setwd('~/git/vespucci-analysis/')
source('R/theme.R')
source('R/functions/utils.R')

# check different resolution
ves_files = c('data/simulations/vespucci/input=circle_res50_50-seed=0-ves_seed=42-max_cells=100.rds','data/simulations/vespucci/input=circle_res200_200-seed=0-ves_seed=42-max_cells=100.rds')
meta_dir = 'data/simulations/objects_meta/'
plot_list = list()
for (ves_file in ves_files) {
    ves_res = readRDS(ves_file)
    meta_file = gsub('-ves_seed.*', '.rds', basename(ves_file))
    input = gsub('input=', '', gsub('-seed=.*', '', meta_file))
    meta_list = readRDS(paste0(meta_dir, meta_file))
    
    truth_pal = pals::kovesi.linear_blue_95_50_c20(100)
    brks = c(1, 2)
    labels = c('min', 'max')
    meta = meta_list$meta
    p1_1 = meta %>%
        mutate(
            label_sim = ifelse(label_sim == 'label1', 'background', label_sim),
            label_sim = factor(label_sim, levels=c('background', 'label2')),
            label_sim = as.integer(label_sim)
        ) %>%
        # arrange(label_sim) %>%
        ggplot(aes(x = y, y = x, fill = label_sim, color = label_sim)) +
        ggtitle(ifelse(grepl('res50_50', ves_file), 'Low resolution', 'High resolution')) +
        ggrastr::rasterise(
            geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_color_gradientn(colours = truth_pal, name = 'Perturbation   ',
                              breaks = brks, labels = labels) +
        scale_fill_gradientn(colours = truth_pal, name = 'Perturbation   ', 
                             breaks = brks, labels = labels) +
        guides(
            fill = guide_colorbar(ticks.colour = NA, frame.colour = NA, title.position='top'),
            colour = guide_colorbar(ticks.colour = NA, frame.colour = NA, title.position='top')
        ) +
        coord_fixed() +
        boxed_theme(size_sm = 5, size_lg = 5) +
        # ggtitle('Perturbation Effect') +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.ticks.length.x = unit(0, 'lines'),
              axis.ticks.length.y = unit(0, 'lines'),
              # legend.position = 'top',
              legend.position = 'none',
              legend.justification = 'bottom',
              legend.key.width = unit(0.18, 'lines'),
              legend.key.height = unit(0.18, 'lines')
        )
    p1_1
    
    aucs = ves_res$spatial_auc_result$aucs
    dat0 = meta %>% left_join(aucs)
    fit = loess(auc ~ x * y, data = dat0, span = 0.015)
    dat0$auc_fit = predict(fit, dat0)
    
    range = range(dat0$auc_fit)
    brks = c(range[1] + 0.1 * diff(range),
             range[2] - 0.1 * diff(range))
    labels = format(range, digits = 2)
    labels = c(paste0(labels[1], ' '),
               paste0(' ', labels[2]))
    
    auc_pal = cet_pal(100, name = 'l19') %>% rev()
    
    p1_2 = dat0 %>%
        # arrange(auc_fit) %>%
        ggplot(aes(x = y, y = x, fill = auc_fit, color = auc_fit)) +
        ggrastr::rasterise(
            geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_fill_gradientn(colours = auc_pal,
                             name = 'AUC   ', labels = labels,
                             limits = range, breaks = brks) +
        scale_color_gradientn(colours = auc_pal,
                              name = 'AUC   ', labels = labels,
                              limits = range, breaks = brks) +
        guides(
            fill = guide_colorbar(ticks.colour = NA, frame.colour = NA, title.position='top'),
            colour = guide_colorbar(ticks.colour = NA, frame.colour = NA, title.position='top')
        ) +
        coord_fixed() +
        boxed_theme(size_lg = 6, size_sm = 5) +
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
              # legend.justification = 'bottom',
              legend.key.width = unit(0.18, 'lines'),
              legend.key.height = unit(0.18, 'lines'),
              plot.title = element_text(size = 6))
    p1_2
    plot_list[[length(plot_list)+1]] = plot_grid(p1_1, p1_2, nrow=2)
}

ves_de_res = readRDS('data/simulations/summaries/de_results/all_vespucci_de_auroc_summary.rds') %>% 
    filter(input %in% c('circle', 'circle_res50_50', 'circle_res200_200'), max_cells == 100, ves_seed == 42) %>% 
    dplyr::select(input, seed, auprc_integral) %>% 
    mutate(de_method = 'vespucci', resolution = factor(case_when(
        input == 'circle_res50_50' ~ 'Low',
        input == 'circle' ~ 'Default',
        input == 'circle_res200_200' ~ 'High'
    ), levels=c('Low', 'Default', 'High')), val = auprc_integral)
labs = ves_de_res %>%
    group_by(resolution) %>%
    summarize(
        # stats_val = median(val),
        stats_val = mean(val),
        val = max(val)
    ) %>%
    ungroup() %>%
    mutate(
        text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 2)),
        text_y = ifelse(val < 0, 0, val)
    )
p2 = ves_de_res %>%
    ggplot(aes(x = resolution, y = val)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA, fill='grey30') +
    geom_label(dat = labs,
               aes(label = text_val, y = text_y), 
               color = ifelse(labs$text_val == 'OOT', 'grey', 'black'),
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0.5, vjust = 0,
               show.legend = FALSE) +
    scale_y_continuous('AUPR', breaks = pretty_breaks(), limits = c(0.18, 1), expand = expansion(c(0.03, 0.125))) +
    scale_x_discrete('Resolution') +
    boxed_theme() +
    theme(aspect.ratio = 1.5,
          axis.text.x = element_text(angle=45, hjust=1),
          legend.position = 'none',
          legend.justification = 'bottom',
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.15, 'lines'))

out_p = plot_grid(wrap_plots(plot_list, nrow=1), p2, rel_widths = c(2,1))
ggsave('fig/final/EFig13/different_resolutions_setup.pdf', out_p, width=8, height=6, units='cm')

spatial_acc_res = readRDS('data/simulations/summaries/spatial_res/spatial_acc.rds') %>% 
    filter((input == 'circle' | grepl('circle_rad', input) | grepl('circle_thick', input)), type == 'AUPRC') %>% 
    dplyr::select(input, seed, val)

spatial_acc_res %<>% 
    mutate(
        input_type = gsub('circle_', '', input),
        rad = as.integer(gsub('_.*', '', gsub('rad', '', input_type))),
        rad = ifelse(is.na(rad), 100, rad),
        thick = as.integer(gsub('.*_', '', gsub('thick', '', input_type))),
        thick = ifelse(is.na(thick), 40, thick)
    ) %>% arrange(rad, thick)
spatial_acc_res %<>% filter(rad <= 100 & thick <= 40, !input %in% c('circle_thick40', 'circle_rad50', 'circle_rad100_thick40'))

# load number of barcodes
meta_files = list.files('data/simulations/objects_meta/', full.names=T)
meta_files_df = map_df(meta_files, convert_filename_to_params) %>% filter(is.na(sp_genes)) %>% type_convert() %>% inner_join(spatial_acc_res %>% dplyr::select(input, seed))
perturbed_barcode_counts = map_df(1:nrow(meta_files_df), function(i){
    print(i)
    tmp = readRDS(meta_files_df$ori_filename[i])$meta
    tmp %<>% 
        mutate(
            label_sim = ifelse(label_sim == 'label1', 'background', label_sim),
            label_sim = factor(label_sim, levels=c('background', 'label2')),
            label_sim = as.integer(label_sim)
        )
    meta_files_df[i,] %>% 
        dplyr::select(-ori_filename) %>%
        mutate(
            no_of_perturbed_barcodes = sum(tmp$label_sim == 2),
            pct_of_perturbed_barcodes = mean(tmp$label_sim == 2)
        )
})

dat0 = spatial_acc_res %>% left_join(perturbed_barcode_counts)

ves_de_res = readRDS('data/simulations/summaries/de_results/all_vespucci_de_auroc_summary.rds') %>% 
    filter((input == 'circle' | grepl('circle_rad', input) | grepl('circle_thick', input)), ves_seed==42, max_cells==100) %>% 
    dplyr::select(input, seed, auprc_integral) 
ves_de_res %<>% 
    mutate(
        input_type = gsub('circle_', '', input),
        rad = as.integer(gsub('_.*', '', gsub('rad', '', input_type))),
        rad = ifelse(is.na(rad), 100, rad),
        thick = as.integer(gsub('.*_', '', gsub('thick', '', input_type))),
        thick = ifelse(is.na(thick), 40, thick)
    ) %>% arrange(rad, thick)
ves_de_res %<>% filter(rad <= 100 & thick <= 40, !input %in% c('circle_thick40', 'circle_rad50', 'circle_rad100_thick40'))

dat1 = ves_de_res %>% left_join(perturbed_barcode_counts)

p3_1 = dat0 %>% 
    ggplot(aes(x=pct_of_perturbed_barcodes, y=val)) +
    # ggplot(aes(x=no_of_perturbed_barcodes, y=val)) +
    ggrastr::rasterise(
        geom_point(size = 0.1, shape='.'), dpi = 600
    ) +
    boxed_theme() +
    scale_x_continuous('% of perturbed barcodes') +
    # scale_x_continuous('# of perturbed barcodes') +
    scale_y_continuous('Spatial AUPRC') +
    theme(aspect.ratio=1)
    # stat_cor(size=1.5, aes(label = ..r.label..))
p3_1

p3_2 = dat1 %>% 
    ggplot(aes(x=pct_of_perturbed_barcodes, y=auprc_integral)) +
    # ggplot(aes(x=no_of_perturbed_barcodes, y=auprc_integral)) +
    ggrastr::rasterise(
        geom_point(size = 0.1, shape='.'), dpi = 600
    ) +
    boxed_theme() +
    scale_x_continuous('% of perturbed barcodes') +
    # scale_x_continuous('# of perturbed barcodes') +
    scale_y_continuous('DE gene AUPRC') +
    theme(aspect.ratio=1)
    # stat_cor(size=1.5, aes(label = ..r.label..))
p3_2

p3 = wrap_plots(p3_1, p3_2, nrow=1)
# ggsave('fig/final/EFig13/no_of_barcodes_against_acc.pdf', p3, width=6, height=5, units='cm')
ggsave('fig/EFig5/pct_of_barcodes_against_acc.pdf', p3, width=6, height=5, units='cm')

perturbed_barcode_counts %>% filter(seed==0, input %in% c('circle', 'circle_rad15_thick25', 'circle_rad50_thick25', 'circle_rad100_thick5'))

ves_files = list.files('data/simulations/vespucci/', full.names=T, pattern='.rds')
ves_files = ves_files[grepl('seed=0', ves_files) & (grepl('rad', ves_files) | grepl('thick', ves_files) | grepl('circle-', ves_files)) & grepl('ves_seed=42-max_cells=100.rds', ves_files)]
ves_files_df = map_df(ves_files, function(ves_file){data.frame(input_file=ves_file, input_type=gsub('input=', '', gsub('.*circle_', '', gsub('-seed.*', '', basename(ves_file)))))}) %>% mutate(
    rad = as.integer(gsub('_.*', '', gsub('rad', '', input_type))),
    rad = ifelse(is.na(rad), 100, rad),
    thick = as.integer(gsub('.*_', '', gsub('thick', '', input_type))),
    thick = ifelse(is.na(thick), 40, thick)
) %>% arrange(rad, thick)
# ves_files_df %<>% filter(input_type %in% c('circle', 'rad5_thick5', 'rad15_thick25', 'rad25_thick40', 'rad50_thick5', 'rad50_thick25', 'rad50_thick40',  'rad100_thick25'))
ves_files_df %<>% filter(input_type %in% c('circle', 'rad15_thick25', 'rad50_thick25', 'rad100_thick5'))
meta_dir = 'data/simulations/objects_meta/'
plot_list = list()
for (i in 1:nrow(ves_files_df)) {
    tmp_row = ves_files_df[i,]
    ves_file = tmp_row$input_file
    ves_res = readRDS(ves_file)
    meta_file = gsub('-ves_seed.*', '.rds', basename(ves_file))
    input = gsub('input=', '', gsub('-seed=.*', '', meta_file))
    meta_list = readRDS(paste0(meta_dir, meta_file))
    truth_pal = pals::kovesi.linear_blue_95_50_c20(100)
    brks = c(1, 2)
    labels = c('min', 'max')
    meta = meta_list$meta
    plot_title = paste0(ifelse(input == 'circle', 'Default\n', ''), paste0('Radius: ', tmp_row$rad, '\nThickness: ', tmp_row$thick))

    meta %<>% 
        mutate(
            label_sim = ifelse(label_sim == 'label1', 'background', label_sim),
            label_sim = factor(label_sim, levels=c('background', 'label2')),
            label_sim = as.integer(label_sim)
        )
    print(table(meta$label_sim))

    p3_1 = meta %>%
        arrange(label_sim) %>%
        ggplot(aes(x = y, y = x, fill = label_sim, color = label_sim)) +
        ggtitle(plot_title) +
        ggrastr::rasterise(
            geom_point(size = 0.3, shape = 21, stroke = 0), dpi = 300
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_color_gradientn(colours = truth_pal, name = 'Perturbation   ',breaks = brks, labels = labels) +
        scale_fill_gradientn(colours = truth_pal, name = 'Perturbation   ', breaks = brks, labels = labels) +
        coord_fixed() +
        boxed_theme(size_sm = 5, size_lg = 5) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length.x = unit(0, 'lines'), axis.ticks.length.y = unit(0, 'lines'), legend.position = 'none') + ggh4x::force_panelsizes(rows = unit(1.1, 'cm'), cols = unit(1.6, 'cm'))
    
    aucs = ves_res$spatial_auc_result$aucs
    dat0 = meta %>% left_join(aucs)
    fit = loess(auc ~ x * y, data = dat0, span = 0.015)
    dat0$auc_fit = predict(fit, dat0)
    
    range = range(dat0$auc_fit)
    brks = c(range[1] + 0.1 * diff(range), range[2] - 0.1 * diff(range))
    labels = format(range, digits = 2)
    labels = c(paste0(labels[1], ' '), paste0(' ', labels[2]))    
    auc_pal = cet_pal(100, name = 'l19') %>% rev()
    
    p3_2 = dat0 %>%
        # arrange(auc_fit) %>%
        ggplot(aes(x = y, y = x, fill = auc_fit, color = auc_fit)) +
        ggrastr::rasterise(geom_point(size = 0.3, shape = 21, stroke = 0), dpi = 300) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_fill_gradientn(colours = auc_pal, name = 'AUC   ', labels = labels,limits = range, breaks = brks) +
        scale_color_gradientn(colours = auc_pal, name = 'AUC   ', labels = labels, limits = range, breaks = brks) +
        guides(fill = guide_colorbar(theme = theme(legend.ticks = element_blank(),legend.frame = element_rect(colour = "black", size = 0.25))), color = guide_colorbar(theme = theme(legend.ticks = element_blank(),legend.frame = element_rect(colour = "black", size = 0.25)))) +
        # ggtitle(tmp_row$input_type) +
        coord_fixed() +
        boxed_theme(size_lg = 6, size_sm = 5) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length.x = unit(0, 'lines'), axis.ticks.length.y = unit(0, 'lines'), legend.position = 'bottom', legend.key.width = unit(0.18, 'lines'), legend.key.height = unit(0.18, 'lines'), plot.title = element_text(size = 6)) + ggh4x::force_panelsizes(rows = unit(1.1, 'cm'), cols = unit(1.6, 'cm'))
    p3_2
    plot_list[[length(plot_list)+1]] = p3_1
    plot_list[[length(plot_list)+1]] = p3_2
}
out_p = wrap_plots(plot_list, nrow=2)
ggsave('fig/EFig5/different_circles_setup.pdf', out_p, width=10, height=8, units='cm')
