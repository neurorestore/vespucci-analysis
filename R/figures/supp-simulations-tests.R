library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(cetcolor)
library(pROC)
library(PRROC)

setwd('~/git/vespucci-analysis/')
source('R/theme.R')

# test different seeds
dat = readRDS('data/simulations/summaries/simulations-auc-summary.rds') %>%
    mutate(
        cor = pearson_cor,
        input_clean = case_when(
            input == 'circle' ~ 'Simulation 1',
            input == 'circle_overlap' ~ 'Simulation 2',
            input == 'stripes' ~ 'Simulation 3',
            input == 'flag' ~ 'Simulation 4'
        )
    )
labs = dat %>%
    group_by(input_clean) %>%
    summarize(
        stats_val = mean(cor),
        val = max(cor)
    ) %>% 
    ungroup() %>%
    mutate(
        text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 2)),
        text_y = ifelse(val < 0, 0, val)
    )
pal = nr_base_4 %>% setNames(unique(dat$input_clean))

p8 = dat %>%
    ggplot(aes(x = input_clean, y = cor)) +
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
p8


dat = readRDS('data/simulations/summaries/simulations-pvals-summary.rds') %>%
    mutate(
        cor = spearman_cor,
        input_clean = case_when(
            input == 'circle' ~ 'Simulation 1',
            input == 'circle_overlap' ~ 'Simulation 2',
            input == 'stripes' ~ 'Simulation 3',
            input == 'flag' ~ 'Simulation 4'
        )
    )
labs = dat %>%
    group_by(input_clean) %>%
    summarize(
        stats_val = mean(cor),
        val = max(cor)
    ) %>% 
    ungroup() %>%
    mutate(
        text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 2)),
        text_y = ifelse(val < 0, 0, val)
    )
pal = nr_base_4 %>% setNames(unique(dat$input_clean))

p9 = dat %>%
    ggplot(aes(x = input_clean, y = cor)) +
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
p9
out_p = wrap_plots(p8, p9, nrow=2)
ggsave('fig/EFig6/simulations-different_seeds.pdf', out_p, width=4, height=6, units='cm')

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
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE,
                                     title.position = "left"),
               color = guide_colorbar(frame.colour = 'black', ticks = FALSE, 
                                      title.position = "left")) +
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
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE),
               color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
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
ggsave('fig/EFig6/different_resolutions_setup.pdf', out_p, width=8, height=6, units='cm')

ves_files = list.files('data/simulations/vespucci/', full.names=T)
ves_files = ves_files[grepl('seed=0', ves_files) & (grepl('rad', ves_files) | grepl('thick', ves_files))]
ves_files_df = map_df(ves_files, function(ves_file){data.frame(input_file=ves_file, input_type=gsub('.*circle_', '', gsub('-seed.*', '', basename(ves_file))))})
ves_files_df %<>% mutate(input_type = factor(input_type, levels=c('rad125', 'rad150', 'thick60', 'thick80', 'rad125_thick60', 'rad125_thick80', 'rad150_thick60', 'rad150_thick80')), perm_name = paste0('Permutation #', as.integer(input_type))) %>% arrange(input_type)
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
    p3_1 = meta %>%
        mutate(
            label_sim = ifelse(label_sim == 'label1', 'background', label_sim),
            label_sim = factor(label_sim, levels=c('background', 'label2')),
            label_sim = as.integer(label_sim)
        ) %>%
        # arrange(label_sim) %>%
        ggplot(aes(x = y, y = x, fill = label_sim, color = label_sim)) +
        ggtitle(tmp_row$perm_name) +
        ggrastr::rasterise(
            geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_color_gradientn(colours = truth_pal, name = 'Perturbation   ',
                              breaks = brks, labels = labels) +
        scale_fill_gradientn(colours = truth_pal, name = 'Perturbation   ', 
                             breaks = brks, labels = labels) +
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE,
                                     title.position = "left"),
               color = guide_colorbar(frame.colour = 'black', ticks = FALSE, 
                                      title.position = "left")) +
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
    
    p3_2 = dat0 %>%
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
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE),
               color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
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
    plot_list[[length(plot_list)+1]] = plot_grid(p3_1, p3_2, nrow=2)
}

out_p = plot_grid(wrap_plots(plot_list, nrow=2))
ggsave('fig/EFig6/different_circles_setup.pdf', out_p, width=8.5, height=9.5, units='cm')


ves_de_res = readRDS('data/simulations/summaries/de_results/all_vespucci_de_auroc_summary.rds') %>% 
    filter((grepl('circle_rad', input) | grepl('circle_thick', input)), max_cells == 100, ves_seed == 42) %>% 
    dplyr::select(input, seed, auprc_integral) %>% 
    mutate(
        de_method = 'vespucci', 
        input_type = factor(gsub('circle_', '', input), levels=c('rad125', 'rad150', 'thick60', 'thick80', 'rad125_thick60', 'rad125_thick80', 'rad150_thick60', 'rad150_thick80')),
        val = auprc_integral,
        perm_name = paste0('Permutation #', as.integer(input_type))
    )
labs = ves_de_res %>%
    group_by(perm_name) %>%
    summarize(
        # stats_val = median(val),
        stats_val = mean(val),
        val = max(val)
    ) %>%
    ungroup() %>%
    mutate(
        text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 3)),
        text_y = ifelse(val < 0, 0, val)
    )
p4 = ves_de_res %>%
    ggplot(aes(x = rev(perm_name), y = val)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA, fill='grey30') +
    geom_label(dat = labs,
               aes(label = text_val, y = text_y), color = 'black',
               # label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0, vjust = 0.5,
               show.legend = FALSE) +
    scale_y_continuous('AUPR', breaks = seq(0.7,1,0.1), limits = c(0.7, 1), expand = expansion(c(0.03, 0.125))) +
    scale_x_discrete('Resolution') +
    boxed_theme() +
    coord_flip() +
    theme(
        # aspect.ratio = 1.8,
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        # axis.text.x = element_text(angle=45, hjust=1),
        legend.position = 'none',
        legend.justification = 'bottom',
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.15, 'lines')
    )
p4
ggsave('fig/EFig6/different_circles_auc.pdf', p4, width=3.7, height=4, units='cm')

# test registration
meta_dir = 'data/simulations/registration/meta/'
meta_files = list.files(meta_dir, full.names=T)
ves_dir = 'data/simulations/registration/vespucci/'

grid_to_plot = data.frame(
    meta_file = meta_files
) %>%
    mutate(
        ves_file = paste0(ves_dir, gsub('\\.rds', '-ves_seed=42-max_cells=100.rds', basename(meta_file))),
        detail = gsub('circle-seed=0-shift=', '', basename(meta_file)),
        shift = gsub('-nreps=.*', '', detail),
        shift_dist = gsub('both', '', shift),
        shift_dist = gsub('right', '', shift_dist),
        shift_dist = gsub('left', '', shift_dist),
        shift_dir = case_when(
            startsWith(shift, 'both') ~ 'both',
            startsWith(shift, 'left') ~ 'left',
            startsWith(shift, 'right') ~ 'right'
        ),
        nreps = gsub('\\.rds',  '', gsub('.*-nreps=', '', detail))
    ) %>%
    dplyr::select(-detail) %>%
    filter(file.exists(ves_file)) %>%
    type_convert() %>%
    arrange(shift_dir, shift_dist, nreps)

grid_to_plot %<>% 
    filter(
        !(!grepl('both', shift) & nreps==3),
        !(grepl('20', shift))
    )

grid_to_plot %<>% mutate(shift_dir = factor(shift_dir, levels=c('left', 'right', 'both')))
grid_to_plot %<>% arrange(shift_dir, shift_dist, nreps) %>% mutate(plot_title = paste0('Permutation #', 8+row_number()))

plot_list = list()
for (i in 1:nrow(grid_to_plot)) {
    tmp_row = grid_to_plot[i,]
    ves_res = readRDS(tmp_row$ves_file)
    meta = readRDS(tmp_row$meta_file)
    
    truth_pal = pals::kovesi.linear_blue_95_50_c20(100)
    brks = c(1, 2)
    labels = c('min', 'max')
    p4_1 = meta %>%
        mutate(
            label_sim = ifelse(label_sim == 'label1', 'background', label_sim),
            label_sim = factor(label_sim, levels=c('background', 'label2')),
            label_sim = as.integer(label_sim),
            alpha = ifelse(label_sim == 1, 0.4, 1)
        ) %>%
        # arrange(label_sim) %>%
        ggplot(aes(x = y, y = x, fill = label_sim, color = label_sim, alpha=alpha)) +
        ggrastr::rasterise(
            geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_color_gradientn(colours = truth_pal, name = 'Perturbation   ',
                              breaks = brks, labels = labels) +
        scale_fill_gradientn(colours = truth_pal, name = 'Perturbation   ', 
                             breaks = brks, labels = labels) +
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE,
                                     title.position = "left"),
               color = guide_colorbar(frame.colour = 'black', ticks = FALSE, 
                                      title.position = "left")) +
        coord_fixed() +
        boxed_theme(size_sm = 5, size_lg = 5) +
        ggtitle(tmp_row$plot_title) +
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
    p4_1
    
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
    
    p4_2 = dat0 %>%
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
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE),
               color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
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
    p4_2
    plot_list[[length(plot_list)+1]] = plot_grid(p4_1, p4_2, nrow=2, align='v')
}

out_p = wrap_plots(plot_list, nrow=2)
ggsave('fig/EFig6/different_registration_setup.pdf', out_p, width=12, height=11, units='cm')

ves_de_res = readRDS('data/simulations/registration/vespucci_de_auroc_summary.rds') %>%
    mutate(
        shift_dist = gsub('both', '', shift),
        shift_dist = gsub('right', '', shift_dist),
        shift_dist = gsub('left', '', shift_dist),
        shift_dir = case_when(
            startsWith(shift, 'both') ~ 'both',
            startsWith(shift, 'left') ~ 'left',
            startsWith(shift, 'right') ~ 'right'
        )
    ) %>%
    arrange(shift_dir, shift_dist, nreps) %>%
    type_convert()

ves_de_res %<>% 
    filter(
        !(!grepl('both', shift) & nreps==3),
        !(grepl('20', shift))
    )

ves_de_res %<>% mutate(shift_dir = factor(shift_dir, levels=c('left', 'right', 'both')))
ves_de_res %<>% arrange(shift_dir, shift_dist, nreps) %>% left_join(grid_to_plot %>% dplyr::select(shift, nreps, plot_title)) %>% mutate(plot_title= factor(plot_title, levels=grid_to_plot$plot_title))
ves_de_res %<>% mutate(val = auprc_integral)
labs = ves_de_res %>%
    group_by(plot_title) %>%
    summarize(
        # stats_val = median(val),
        stats_val = mean(val),
        val = max(val)
    ) %>%
    ungroup() %>%
    mutate(
        text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 3)),
        text_y = ifelse(val < 0, 0, val)
    )
p5 = ves_de_res %>%
    ggplot(aes(x = plot_title, y = val)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA, fill='grey30') +
    geom_label(dat = labs,
               aes(label = text_val, y = text_y), color = 'black',
               # label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0, vjust = 0.5,
               show.legend = FALSE) +
    scale_y_continuous('AUPR', breaks = seq(0.7,1,0.1), limits = c(0.7, 1), expand = expansion(c(0.03, 0.125))) +
    scale_x_discrete('Resolution') +
    boxed_theme() +
    coord_flip() +
    theme(
        # aspect.ratio = 1.8,
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        # axis.text.x = element_text(angle=45, hjust=1),
        legend.position = 'none',
        legend.justification = 'bottom',
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.15, 'lines')
    )
p5
ggsave('fig/EFig6/different_registration_auc.pdf', p5, width=4, height=10, units='cm')


# plot different subsampling
meta = readRDS('data/simulations/objects_meta/input=circle-seed=0.rds')$meta
max_cells = c(10, 50, 100, 200, 500, 1000)
plot_list = list()
for (curr_max_cells in max_cells) {
    ves_res = readRDS(paste0('data/simulations/vespucci/input=circle-seed=0-ves_seed=42-max_cells=',curr_max_cells,'.rds'))
    
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
    
    p6_1 = dat0 %>%
        # arrange(auc_fit) %>%
        ggplot(aes(x = y, y = x, fill = auc_fit, color = auc_fit)) +
        ggrastr::rasterise(
            geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
        ) +
        ggtitle(paste0('Max barcodes: ', curr_max_cells)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_fill_gradientn(colours = auc_pal,
                             name = 'AUC   ', labels = labels,
                             limits = range, breaks = brks) +
        scale_color_gradientn(colours = auc_pal,
                              name = 'AUC   ', labels = labels,
                              limits = range, breaks = brks) +
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE),
               color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
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
    p6_1
    plot_list[[length(plot_list)+1]] = p6_1
}
out_p = wrap_plots(plot_list, nrow=2)
ggsave('fig/EFig6/different_subsample_setup.pdf', out_p, width=6.5, height=6.5, units='cm')


ves_de_res = readRDS('data/simulations/summaries/de_results/all_vespucci_de_auroc_summary.rds') %>% 
    filter(input == 'circle', ves_seed==42) %>% 
    dplyr::select(input, seed, auprc_integral, max_cells) %>% 
    mutate(plot_title = paste0('Max barcodes: ', max_cells), val=auprc_integral) %>% arrange(max_cells)
ves_de_res %<>% mutate(plot_title=factor(plot_title, levels=unique(ves_de_res$plot_title)))
labs = ves_de_res %>%
    group_by(plot_title) %>%
    summarize(
        # stats_val = median(val),
        stats_val = mean(val),
        val = max(val)
    ) %>%
    ungroup() %>%
    mutate(
        text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 3)),
        text_y = ifelse(val < 0, 0, val)
    )
p7 = ves_de_res %>%
    ggplot(aes(x = plot_title, y = val)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA, fill='grey30') +
    geom_label(dat = labs,
               aes(label = text_val, y = text_y), color = 'black',
               # label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0, vjust = 0.5,
               show.legend = FALSE) +
    scale_y_continuous('AUPR', breaks = seq(0.7,1,0.1), limits = c(0.7, 1), expand = expansion(c(0.03, 0.125))) +
    scale_x_discrete('Resolution') +
    boxed_theme() +
    coord_flip() +
    theme(
        # aspect.ratio = 1.8,
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        # axis.text.x = element_text(angle=45, hjust=1),
        legend.position = 'none',
        legend.justification = 'bottom',
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.15, 'lines')
    )
p7
ggsave('fig/EFig6/different_subsample_auc.pdf', p7, width=4, height=4,5, units='cm')
