setwd('~/git/vespucci-analysis')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(lawstat)
library(nparcomp)
library(ontologyIndex)
source('R/theme.R')


###############################################################################-
## a. RCTD ####
###############################################################################-

meta = readRDS('data/published_data/seurat/Calcagno2022.rds')@meta.data %>%
    mutate(
        barcode = gsub('-', '_', barcode),
        label = timepoint
    )
rctd_res = readRDS('data/published_data/rctd/Calcagno2022/rctd-summary-res.rds') %>%
    mutate(barcode = gsub('-', '_', barcode))

meta %<>% left_join(
    rctd_res,
    by='barcode'    
) %>%
    filter(!is.na(cell_type)) %>%
    dplyr::rename(cell_type_marker = cell_type)

cell_type_markers = unique(meta$cell_type_marker)
cell_type_mapping = data.frame(
    cell_type_marker = c(
        c('Gsn', 'Postn', 'Cxcl5'), # Cxcl5 from paper
        'Myo1b',
        'Myh11',
        'Rep',
        c('DC', 'Mono', 'Arg1', 'NSG', 'Lyve1'),
        c('SigF', 'Retnlg', 'ISG'),
        c('EC', 'Npr3', 'Fbln5', 'Kit', 'Pecam1'), # Kit, Pecam1 from EC
        c('Myh6', 'Ankrd1', 'Xirp2')
    ),
    sub_cell_type = c(
        rep('Fibroblasts', 3),
        'Smooth muscle cells',
        # 'Myh11',
        # 'Rep',
        c('Others', 'Others'),
        rep('Macrophages', 5),
        rep('Neutrophils', 3),
        rep('Endothelial cells', 5),
        rep('Cardiomyoctes', 3)
    )
)

meta %<>%
    left_join(cell_type_mapping)
table(meta$sub_cell_type)
sum(is.na(meta$sub_cell_type))

# change neutrophils (only 1) and NA to others
meta %<>%
    mutate(
        sub_cell_type = ifelse(is.na(sub_cell_type), 'Others', sub_cell_type),
        sub_cell_type = ifelse(sub_cell_type == 'Neutrophils', 'Others', sub_cell_type),
    )
table(meta$sub_cell_type)
sum(is.na(meta$sub_cell_type))

color_pal = c(
    'Cardiomyoctes' = "#A7C2A2", 
    'Others' = "#FFDCA4",
    'Endothelial cells' = "#ADDEDA", 
    'Fibroblasts' = "#887AA1", 
    'Macrophages' = "#F69DA4"
)

p1 = meta %>%
    ggplot(aes(x=ori_x,y=ori_y, fill=sub_cell_type, color=sub_cell_type)) +
    ggrastr::rasterise(
        geom_point(size = 0.4, shape = 21, stroke = 0), dpi = 300
    ) +
    boxed_theme() +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0, 'lines'),
          axis.ticks.length.y = unit(0, 'lines'),
          legend.position = 'none',
          strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    facet_wrap(~replicate, ncol = 6, scales='free')
# p1

p2 = meta %>%
    ggplot(aes(x=x,y=y, fill=sub_cell_type, color=sub_cell_type)) +
    ggrastr::rasterise(
        geom_point(size = 0.4, shape = 21, stroke = 0), dpi = 300
    ) +
    boxed_theme() +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0, 'lines'),
          axis.ticks.length.y = unit(0, 'lines'),
          legend.position = 'right',
          legend.justification = 'bottom',
          legend.text = element_text(size=5),
          legend.key.size = unit(0.4, 'lines'),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=0.8))) +
    facet_wrap(~replicate, ncol = 6, scales='free')
p0 = wrap_plots(p1,p2, nrow=2)
# p0
ggsave('fig/final/EFig4/rctd-registration.pdf', p0, width = 11, height = 4, units='cm')

###############################################################################-
## b. genes ####
###############################################################################-

# read Calcagno dataset
sc = readRDS('data/published_data/seurat/Calcagno2022.rds')
meta = sc@meta.data %>%
  mutate(
    barcode = gsub('-', '_', barcode),
    label = timepoint
  )
expr = GetAssayData(sc, slot = 'counts')
colnames(expr) = gsub('-', '_', colnames(expr))

# normalize
expr %<>% NormalizeData()
# extract coordinates
dat0 = meta %>% dplyr::select(barcode, x, y, label)

# list genes to plot
genes_to_plot = c(
  # up in d7
  'GPNMB', 'PTN', 'ITGBL1',
  # up in d1
  'HOPX', 'HRC',
  'HMOX1', 'SPP1'
)
stopifnot(all(str_to_title(genes_to_plot) %in% rownames(expr)))

# iterate through genes
expr_plots = list()
for (idx in seq_along(genes_to_plot)) {
  gene = genes_to_plot[idx]
  dat0$expr = expr[str_to_title(gene), dat0$barcode]
  
  # extract data frame
  plot_df = data.frame()
  conditions = unique(dat0$label)
  for (condition in conditions) {
    tmp_plot_df = dat0 %>% dplyr::filter(label == !!condition)
    fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.02, degree = 1)
    tmp_plot_df$expr = predict(fit, tmp_plot_df)
    plot_df %<>% rbind(tmp_plot_df)
  }
  
  range = range(plot_df$expr)
  brks = c(range[1] + 0.1 * diff(range),
           range[2] - 0.1 * diff(range))
  labels = c('min', 'max')
  
  pal = scico::scico(palette = 'davos', n = 100)
  pal = jdb_palette("horizon") %>% as.character()
  pal = pals::kovesi.linear_bmy_10_95_c71(100)
  pal = pals::kovesi.linear_blue_5_95_c73(100)
  pal = nr_heat_red_spatial
  pal = pals::kovesi.linear_grey_10_95_c0(100)
  pal = pals::kovesi.diverging_isoluminant_cjm_75_c23(100)
  pal = pals::ocean.dense(100)
  pal = brewer.pal(9, 'RdGy') %>% rev
  pal = pals::ocean.thermal(100)
  pal = nr_heat_red_no_white %>% tail(-5)
  p = plot_df %>%
      mutate(label = ifelse(label == 'd1', 'Day 3', 'Day 7')) %>%
      # arrange(-expr) %>%
      ggplot(aes(x = x, y = y, fill = expr)) +
      rasterise(geom_point(size = 0.3,
                           shape = 21, stroke = 0, alpha = 1), dpi = 600) +
      scale_color_gradientn(colours = pal, name = 'Expression',
                            breaks = brks, labels = labels) +
      scale_fill_gradientn(colours = pal, name = 'Expression', 
                           breaks = brks, labels = labels) +
      guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
      boxed_theme(size_lg = 6, size_sm = 5) +
      ggtitle(toupper(gene)) +
      theme(
          aspect.ratio = 1,
          plot.title = element_text(size = 5, margin = margin(0,0,-5,0)),
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.1), "cm"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.14, 'lines'),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          legend.position = 'right',
          legend.justification = 'bottom'
      ) +
      facet_wrap(~label) +
      coord_fixed()
  if (idx != 10)
    p = p + theme(legend.position = 'none')
  p
  expr_plots[[idx]] = p
}
p1 = wrap_plots(expr_plots, ncol = 5)
ggsave('fig/final/EFig4/genes.pdf', p1, width = 10, height = 4.5,
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## c. GO terms ####
###############################################################################-

sc = readRDS('data/published_data/seurat_GO/Calcagno2022.rds')
go_df = readRDS('data/metadata/go_names.rds')
meta = sc@meta.data %>%
    mutate(barcode = gsub('-', '_', barcode))
mat = GetAssayData(sc, slot='counts')
gos_to_plot = c(
    'age-dependent response to reactive oxygen species',
    'myeloid leukocyte activation',
    'NADH dehydrogenase activity'
)
gos = go_df %>% filter(go_name %in% gos_to_plot) %>%
    pull(go) %>%
    gsub('\\:', '-', .)
mat = mat[gos,]

# iterate through GO terms
go_plots = expr_plots
for (idx in seq_len(nrow(mat))) {
  go = rownames(mat)[idx]
  title = go_df$go_name[gsub('\\:', '-', go_df$go) == go]
    
  # extract data frame
  plot_df = data.frame()
  conditions = unique(meta$label)
  meta$expr = mat[go, meta$barcode]
  for (condition in conditions) {
    tmp_plot_df = meta %>% 
      filter(label == !!condition)
    fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.02, degree = 1)
    tmp_plot_df$expr = predict(fit, tmp_plot_df)
    plot_df %<>% rbind(tmp_plot_df)
  }
  
  range = range(plot_df$expr)
  brks = c(range[1] + 0.1 * diff(range),
           range[2] - 0.1 * diff(range))
  labels = c('min', 'max')
  
  pal = pals::ocean.thermal(100)
  pal = brewer.pal(9, 'RdGy') %>% rev
  pal = pals::ocean.solar(100)
  # pal = nr_heat_blue_no_white %>% tail(-5)
  pal = nr_heat_blue_spatial
  
  title = paste(strwrap(title, 30), collapse = '\n')
  
  p = plot_df %>%
      mutate(label = ifelse(label == 'd1', 'Day 3', 'Day 7')) %>%
      # arrange(-expr) %>%
      ggplot(aes(x = x, y = y, fill = expr)) +
      rasterise(geom_point(size = 0.3,
                           shape = 21, stroke = 0, alpha = 1), dpi = 600) +
      scale_color_gradientn(colours = pal, name = 'Expression',
                            breaks = brks, labels = labels) +
      scale_fill_gradientn(colours = pal, name = 'Expression', 
                           breaks = brks, labels = labels) +
      guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
      boxed_theme(size_lg = 6, size_sm = 5) +
      ggtitle(title) +
      theme(
          aspect.ratio = 1,
          plot.title = element_text(size = 5, margin = margin(0,0,-5,0)),
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.1), "cm"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.14, 'lines'),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          legend.position = 'right',
          legend.justification = 'bottom'
      ) +
      facet_wrap(~label) +
      coord_fixed()
  if (idx != 10)
    p = p + theme(legend.position = 'none')
  p
  go_plots[[length(go_plots)+1]] = p
}
p2 = wrap_plots(go_plots, nrow = 2)
ggsave('fig/final/EFig4/genes-GO-modules.pdf', p2, width = 16, height = 5,
       units = 'cm', useDingbats = FALSE)


###############################################################################-
## d. lollipop plot: genes ####
###############################################################################-

# load results
dat0 = readRDS('data/published_data/vespucci/Calcagno2022-seed=42-nsub=10.rds')$de_feature_result
dat0$p_val = vapply(dat0$p_val, function(x){max(min_pval, x)}, as.numeric(1))
dat0 %<>%
  mutate(up_reg = effect_size > 0) %>%
  group_by(up_reg) %>%
  arrange(p_val, -abs(effect_size)) %>%
  mutate(
    log_p_val = -log10(p_val)
  ) %>%
  slice(1:15)
dat0 %<>% 
  arrange(-effect_size)
dat0$gene = factor(dat0$gene, levels=rev(as.character(dat0$gene)))

p3 = dat0 %>%
  ggplot(aes(x = gene, y = effect_size)) +
  # facet_wrap(sign(effect_size) ~ ., ncol = 1, scales = 'free') +
  geom_hline(aes(yintercept = 0), size = 0.4, color = 'grey88') +
  geom_segment(aes(xend = gene, yend = 0), color = 'grey88') +
  geom_point(shape = 21, stroke = 0.2, size = 0.9, color = 'black', 
             fill = 'grey80') +
  geom_segment(aes(y = -Inf, yend = Inf, x = 15.5, xend = 15.5), 
               color = 'grey20', size = 0.225) +
  scale_y_continuous(expression(beta)) +
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
p3
ggsave('fig/final/EFig4/lollipop-genes.pdf', p3, width = 8, height = 8, 
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## e. lollipop plot: GO modules ####
###############################################################################-

# load GO
dat0 = readRDS('data/published_data/vespucci_GO/Calcagno2022-seed=42-nsub=10.rds')$de_feature_result
dat0$p_val = vapply(dat0$p_val, function(x){max(min_pval, x)}, as.numeric(1))
dat0 %<>%
  mutate(up_reg = effect_size > 0) %>%
  group_by(up_reg) %>%
  arrange(p_val, -abs(effect_size)) %>%
  mutate(
    log_p_val = -log10(p_val)
  ) %>%
  slice(1:15)

# name GO terms
dat0 %<>%
  mutate(gene = chartr('-', ':', gene)) %>% 
  mutate(descr = go$name[gene]) 
dat0 %<>% 
  arrange(-effect_size)
dat0$descr = factor(dat0$descr, levels=rev(as.character(dat0$descr)))

# plot
p4 = dat0 %>%
  ggplot(aes(x = descr, y = effect_size)) +
  # facet_wrap(sign(effect_size) ~ ., ncol = 1, scales = 'free') +
  geom_hline(aes(yintercept = 0), size = 0.4, color = 'grey88') +
  geom_segment(aes(xend = descr, yend = 0), color = 'grey88') +
  geom_point(shape = 21, stroke = 0.2, size = 0.9, color = 'black', 
             fill = 'grey80') +
  geom_segment(aes(y = -Inf, yend = Inf, x = 15.5, xend = 15.5), 
               color = 'grey20', size = 0.225) +
  scale_y_continuous(expression(beta)) +
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
p4
ggsave('fig/final/EFig4/lollipop-GO.pdf', p4, width = 18, height = 8, 
       units = 'cm', useDingbats = FALSE)

