library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(cetcolor)
library(pROC)
library(PRROC)

# setwd('C:/Users/teo/Documents/EPFL/projects/vespucci/')
setwd('~/git/vespucci/')
source('R/theme.R')

# test registration
# meta_dir = 'data/rejected_review/test_registration/simulations/meta/'
meta_dir = '/work/upcourtine/vespucci/rejected_review/test_registration/simulations/meta/'
meta_files = list.files(meta_dir, full.names=T)
# ves_dir = 'data/rejected_review/test_registration/simulations/vespucci/'
ves_dir = '/work/upcourtine/vespucci/rejected_review/test_registration/simulations/vespucci/'

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
grid_to_plot %<>% filter(!is.na(shift_dir))
grid_to_plot %<>% arrange(shift_dir, shift_dist, nreps) %>% mutate(plot_title = paste0('Permutation #', row_number()))

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
        ggrastr::rasterise(geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_color_gradientn(colours = truth_pal, name = 'Perturbation   ',breaks = brks, labels = labels) +
        scale_fill_gradientn(colours = truth_pal, name = 'Perturbation   ',breaks = brks, labels = labels) +
        guides(color = guide_colorbar(theme = theme(legend.ticks = element_blank(),legend.frame = element_rect(colour = "black", size = 0.25)))) +
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
    # p4_1
    aucs = ves_res$spatial_auc_result$aucs
    dat0 = meta %>% left_join(aucs)
    fit = loess(auc ~ x * y, data = dat0, span = 0.015)
    dat0$auc_fit = predict(fit, dat0)
    range = range(dat0$auc_fit)
    brks = c(range[1] + 0.1 * diff(range),range[2] - 0.1 * diff(range))
    labels = format(range, digits = 2)
    labels = c(paste0(labels[1], ' '),paste0(' ', labels[2]))
    auc_pal = cet_pal(100, name = 'l19') %>% rev()
    p4_2 = dat0 %>%
        # arrange(auc_fit) %>%
        ggplot(aes(x = y, y = x, fill = auc_fit, color = auc_fit)) +
        ggrastr::rasterise(geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_fill_gradientn(colours = auc_pal,name = 'AUC   ', labels = labels,limits = range, breaks = brks) +
        scale_color_gradientn(colours = auc_pal, name = 'AUC   ', labels = labels,limits = range, breaks = brks) +
        guides(fill = guide_colorbar(theme = theme(legend.ticks = element_blank(),legend.frame = element_rect(colour = "black", size = 0.25)))) +
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
    # p4_2
    plot_list[[length(plot_list)+1]] = plot_grid(p4_1, p4_2, nrow=2, align='v')
}
out_p = wrap_plots(plot_list, nrow=2)
ggsave('fig/final/EFig11/different_registration_setup.pdf', out_p, width=12, height=11, units='cm')

ves_de_res = readRDS('data/rejected_review/test_registration/simulations/vespucci_de_auroc_summary.rds') %>%
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
ggsave('fig/final/EFig11/simdifferent_registration_auc.pdf', p5, width=4, height=6, units='cm')

meta_dir = 'data/rejected_review/test_registration/real_data/meta/'
meta_files = list.files(meta_dir, full.names=T, pattern='Calcagno2022')
ves_dir = 'data/rejected_review/test_registration/real_data/vespucci/'
grid_to_plot = data.frame(
  meta_file = meta_files
) %>%
  mutate(
    ves_file = paste0(ves_dir, gsub('\\.rds', '-seed=42-nsub=10.rds', basename(meta_file))),
    detail = gsub('.*shift=', '', basename(meta_file)),
    shift = gsub('-nreps=.*', '', detail),
    shift_dist = gsub('both', '', shift),
    shift_dist = gsub('right', '', shift_dist),
    shift_dist = gsub('left', '', shift_dist),
    shift_dir = case_when(
      startsWith(shift, 'both') ~ 'both',
      startsWith(shift, 'left') ~ 'left',
      startsWith(shift, 'right') ~ 'right'
    ),
    nreps = gsub('\\.rds',  '', gsub('.*-nreps=', '', detail)),
    label2 = ifelse(grepl('treated_old', ves_file), 'treated', 'young')
  ) %>%
  dplyr::select(-detail) %>%
  filter(file.exists(ves_file)) %>%
  type_convert()

annot_factors = c('Remote zone', 'Border zone 1', 'Border zone 2', 'Infarction zone')
color_pal = c('#746FB4', '#AB4567', '#E3201D', '#BEBEBE') %>% set_names(annot_factors)

grid_to_plot %<>% filter(shift_dir == 'both', nreps==2)
shifts = paste0('both', seq(20, 180, 20))

plot_list = list()
for (shift in shifts){
  print(shift)
  tmp_row = grid_to_plot %>% filter(shift == !!shift)
  
  ves_res = readRDS(tmp_row$ves_file)[[1]]
  aucs = ves_res$spatial_auc_result$aucs
  meta = readRDS(tmp_row$meta_file) %>% mutate(barcode = gsub('-', '_', barcode), annotation_name = factor(annotation_name, levels=annot_factors))
  dat0 = aucs %>% left_join(meta)
  
  p1a = dat0 %>%
    mutate(
        label_clean = case_when(label == 'd1' ~ 'Day 3', T ~ 'Day 7'),
        replicate_clean = gsub('rep','Rep',gsub('d7','Day 7', gsub('d1', 'Day 3', gsub('_', ' ', gsub('MI_', '', replicate))))),
    ) %>%
    filter(label_clean == 'Day 7') %>%
    ggplot(aes(x=x,y=y, fill=annotation_name, color=annotation_name)) +
    ggrastr::rasterise(
      geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 600
    ) +
    ggtitle(paste0('Permutation #', which(shift == shifts))) +
    boxed_theme() +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    theme(
      aspect.ratio=1,
      plot.margin = unit(c(0.05, 0.05, 0.05, 0.5), "cm"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_blank(),
      axis.ticks.length.x = unit(0, 'lines'),
      axis.ticks.length.y = unit(0, 'lines'),
      legend.position = 'none',
      legend.justification = 'bottom',
      legend.text = element_text(size=5),
      legend.title = element_blank(),
      strip.background = element_blank(),
      # strip.text.x = element_blank(),
      legend.spacing.x = unit(0.01, 'cm')
    ) + 
    facet_wrap(~replicate_clean, nrow=1, scales='free') +
    guides(colour = guide_legend(override.aes = list(size=1)))
  
  fit = loess(auc ~ x * y, data = dat0, span = 0.015)
  dat0$auc_fit = predict(fit, dat0)
  
  range = range(dat0$auc_fit)
  brks = c(range[1] + 0.1 * diff(range),
           range[2] - 0.1 * diff(range))
  labels = format(range, digits = 2)
  # labels = c(paste0(labels[1], ' '), paste0(' ', labels[2]))
  
  auc_pal = cet_pal(100, name = 'l19') %>% rev()
  p1b = dat0 %>%
    # arrange(auc_fit) %>%
    ggplot(aes(x = x, y = y, fill = auc_fit, color = auc_fit)) +
    ggrastr::rasterise(
      geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 600
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
    # coord_fixed() +
    boxed_theme(size_lg = 5, size_sm = 5) +
    theme(
      plot.margin = unit(c(0.05, 0.05, 0.05, 0.1), "cm"),
      aspect.ratio=1,
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
      legend.position = 'none',
      legend.justification = 'right',
      legend.key.width = unit(0.15, 'lines'),
      legend.key.height = unit(0.15, 'lines')
      # plot.title = element_blank()
    )
  p1_0 = wrap_plots(p1a, p1b, widths=c(3,1))
  p1_0
  plot_list[[length(plot_list)+1]] = p1_0
}

final_p = wrap_plots(plot_list, nrow=3)
ggsave(paste0('fig/final/EFig11/Calcagno2022-registered-aucs.pdf'), final_p, height=8, width=14, units='cm')

meta_dir = 'data/rejected_review/test_registration/real_data/meta/'
meta_files = list.files(meta_dir, full.names=T, pattern='regen_final')
ves_dir = 'data/rejected_review/test_registration/real_data/vespucci/'
grid_to_plot = data.frame(
  meta_file = meta_files
) %>%
  mutate(
    ves_file = paste0(ves_dir, gsub('\\.rds', '-seed=42-nsub=10.rds', basename(meta_file))),
    detail = gsub('.*shift=', '', basename(meta_file)),
    shift = gsub('-nreps=.*', '', detail),
    shift_dist = gsub('both', '', shift),
    shift_dist = gsub('right', '', shift_dist),
    shift_dist = gsub('left', '', shift_dist),
    shift_dir = case_when(
      startsWith(shift, 'both') ~ 'both',
      startsWith(shift, 'left') ~ 'left',
      startsWith(shift, 'right') ~ 'right'
    ),
    nreps = gsub('\\.rds',  '', gsub('.*-nreps=', '', detail)),
    label2 = ifelse(grepl('treated_old', ves_file), 'treated', 'young')
  ) %>%
  dplyr::select(-detail) %>%
  filter(file.exists(ves_file)) %>%
  type_convert()
fixed_coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5))

grid_to_plot %<>% filter(shift_dir == 'both', nreps==2)
shifts = paste0('both', c(50, 100, 150, 200))

plot_list = list()
for (shift in shifts){
  print(shift)
  grid_to_plot2 = grid_to_plot %>% filter(shift == !!shift)
  meta = map_df(grid_to_plot2$meta_file, function(meta_file) {
    readRDS(meta_file) %>% mutate(barcode=gsub('-', '_', barcode)) 
  }) %>% distinct()
    
  conditions = c('young', 'old', 'treated')
  coords = map_dfr(conditions, ~ {
    condition = .x
    fit = loess(expr ~ x * y, data = filter(meta, label == condition),
                span = 0.02, degree = 1)
    coords = fixed_coords %>%
      mutate(interp = predict(fit, newdata = .)) %>%
      drop_na(interp) %>%
      mutate(condition = condition)
  })
  # range = range(coords$interp)
  range = quantile(coords$interp, c(0.01, 0.99))
  brks = c(range[1] + 0.1 * diff(range),
           range[2] - 0.1 * diff(range))
  labels = c('min', 'max')
  p1a = coords %>%
    mutate(
      interp = winsorize(interp, quantile(coords$interp, c(0.01, 0.99))),
      gene = 'Gfap',
      condition = factor(str_to_title(condition), levels=c('Old', 'Young', 'Treated'))
    ) %>%
    ggplot(aes(x = x, y = y, fill = interp)) +
    facet_grid(gene~ condition, switch='y') +
    geom_raster() +
    scale_y_reverse(expand = c(0, 0), breaks = c(-150, 150),
                    labels = c(expression('R'%->%""),
                               expression(""%<-%'L')),
                    name = 'Gfap'
    ) +
    scale_x_continuous(expand = c(0, 0), breaks = c(-400, 400),
                       labels = c(expression(""%<-%'Rostral'),
                                  expression('Caudal'%->%""))
    ) +
    scale_fill_gradientn(name = 'Expression', colors = nr_tree_red,
                         limits = range, breaks = brks, labels = labels) +
    guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
    coord_fixed() +
    boxed_theme() +
    ggtitle(paste0('Permutation #', which(shift == shifts))) +
    theme(
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
      legend.key.width = unit(0.18, 'lines'),
      legend.key.height = unit(0.18, 'lines')
    )
  p1a
  
  plot_list2 = list()
  for (curr_label in c('young', 'treated')) {
    
    auc_df = readRDS(grid_to_plot2 %>% filter(label2 == curr_label) %>% pull(ves_file))[[1]]$spatial_auc_result$aucs %>% inner_join(meta)
    fit = loess(auc ~ x * y, data = auc_df, span = 0.02, degree = 1)
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
    
    p1b = dat1 %>%
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
      guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE),
             color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
      coord_fixed() +
      boxed_theme() +
      ggtitle(paste0(str_to_title(curr_label), ' vs. old')) +
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
    plot_list2[[length(plot_list2)+1]] = p1b
  }
    
  out_p = wrap_plots(p1a, wrap_plots(plot_list2), nrow=1, widths=c(3,2))
  out_p
  plot_list[[length(plot_list)+1]] = out_p
}

final_p = wrap_plots(plot_list, ncol=1)
ggsave(paste0('fig/final/EFig11/regen_final-registered-aucs.pdf'), final_p, height=12, width=10, units='cm')
