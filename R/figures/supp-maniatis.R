setwd('~/git/vespucci-analysis/')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(lawstat)
library(nparcomp)
library(ontologyIndex)
source('R/theme.R')

###############################################################################-
## a. AUC ####
###############################################################################-

# read dataset
sc = readRDS('data/real_data/seurat/Maniatis2019_subset.rds')
meta = sc@meta.data %>%
  mutate(
    x = ori_x
  )

ves_res0 = readRDS('data/real_data/vespucci/Maniatis2019-seed=42-nsub=10.rds')
ves_de_res0 = readRDS('data/real_data/DE_summaries/vespucci/Maniatis2019-seed=42-nsub=10-de=nebula_nbgmm.rds') %>%
  arrange(comparison, p_val, -abs(logFC))

comparison = 'WT-SOD|p100'
ves_res = ves_res0[[comparison]]$aucs
ves_de_res = ves_de_res0 %>%
  filter(comparison == !!comparison)

dat0 = meta %>%
  left_join(ves_res) %>%
  filter(!is.na(auc))

fit = loess(auc ~ x * y, data = dat0, span = 0.02, degree = 1)
dat0$auc_fit = predict(fit, dat0)

range = range(dat0$auc_fit)
brks = c(range[1] + 0.1 * diff(range),
         range[2] - 0.1 * diff(range))
labels = format(range, digits = 2)

library(cetcolor)
auc_pal = pals::kovesi.linear_kryw_5_100_c67(100) %>% rev %>% tail(-5)
auc_pal = nr_heat_red_spatial
auc_pal = nr_heat_red_no_white
auc_pal = cet_pal(100, name = 'l19') %>% rev()
# auc_pal = cet_pal(100, name = 'l18') # %>% rev()
p1 = dat0 %>%
  # arrange(auc_fit) %>%
  # mutate(auc_fit = winsorize(auc_fit, c(0.55, NA))) %>% 
  ggplot(aes(x = x, y = y, fill = auc_fit, color = auc_fit)) +
  ggrastr::rasterise(
    geom_point(size = 0.4, shape = 21, stroke = 0), dpi = 300
  ) +
  scale_fill_gradientn(colours = auc_pal,
                       name = 'AUC', labels = labels,
                       breaks = brks) +
  scale_color_gradientn(colours = auc_pal,
                        name = 'AUC', labels = labels,
                        breaks = brks) +
  guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE),
         color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
  coord_fixed() +
  boxed_theme() +
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
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.15, 'lines'),
        panel.background = element_rect(fill = 'grey96'),
        plot.title = element_text(size = 5))
p1
ggsave("fig/EFig9/AUC.pdf", p1, width = 5, height = 5, units = 'cm',
       useDingbats = FALSE)

###############################################################################-
## b. genes ####
###############################################################################-

meta = sc@meta.data %>%
  mutate(
    x = ori_x,
    label = ifelse(grepl('WT', label_ori), 'WT', 'SOD'),
    label = factor(label, levels = c('WT', 'SOD'))
  )
normalize_data = TRUE
expr = GetAssayData(sc, slot = 'counts')
colnames(expr) = gsub('-', '_', colnames(expr))

if (normalize_data) {
  expr %<>% NormalizeData()
  norm_suffix = '-norm'
} else {
  norm_suffix = '-raw'    
}

# extract coordinates
dat0 = meta %>% dplyr::select(barcode, x, y, label)

# list genes to plot
genes_to_plot = c(
  'NEFL', 
  'NEFH', 
  'NEFM',
  'SOD1',
  'MT1',
  'CST3'
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
  pal = nr_heat_yellow_spatial
  pal = nr_heat_red_no_white %>% tail(-5)
  p = plot_df %>%
    # mutate(label = ifelse(label == 'd1', 'Day 3', 'Day 7')) %>%
    # arrange(-expr) %>%
    ggplot(aes(x = x, y = y, fill = expr)) +
    rasterise(geom_point(size = 0.4,
                         shape = 21, stroke = 0, alpha = 1), dpi = 600) +
    scale_color_gradientn(colours = pal, name = 'Expression',
                          breaks = brks, labels = labels) +
    scale_fill_gradientn(colours = pal, name = 'Expression', 
                         breaks = brks, labels = labels) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    boxed_theme(size_lg = 6, size_sm = 5) +
    ggtitle(toupper(gene)) +
    coord_fixed() +
    theme(
      # aspect.ratio = 1,
      plot.title = element_text(size = 5, 
                                margin = margin(0,0,-2,0)),
      plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
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
      legend.justification = 'bottom',
      panel.background = element_rect(fill = 'grey96')
    ) +
    facet_wrap(~label)
  if (idx != length(genes_to_plot))
    p = p + theme(legend.position = 'none')
  p
  expr_plots[[idx]] = p
}
p1 = wrap_plots(expr_plots, ncol = 2)
ggsave('fig/EFig9/genes.pdf', p1, width = 10, height = 8,
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## c. GO terms ####
###############################################################################-

# load GO
go = get_ontology("data/GO/go.obo")

# read GO sub-matrix
# mat = readRDS('data/real_data/seurat_GO/Maniatis2019_subset_GO.rds')
mat = readRDS('data/real_data/seurat_GO/DE/Maniatis2019.rds')[[1]]
mat = GetAssayData(mat, slot='count')
meta = data.frame(barcode = colnames(mat))
# merge in coordinates from single-cell object
coords = sc@meta.data %>% 
  dplyr::select(barcode, label_ori, ori_x, y) %>% 
  mutate(x = ori_x,
         label = ifelse(grepl('WT', label_ori), 'WT', 'SOD'),
         label = factor(label, levels = c('WT', 'SOD')))
meta %<>% left_join(coords, by = 'barcode')

go_df = readRDS('data/metadata/go_names.rds')
genes_to_plot = c(
    'NEUROFILAMENT BUNDLE ASSEMBLY',
    'NEURON PROJECTION',
    'PERIPHERAL NERVOUS SYSTEM AXON REGENERATION',
    'G PROTEIN-COUPLED RECEPTOR SIGNALING PATHWAY'
)

# iterate through GO terms
go_plots = list()
for (gene in genes_to_plot) {
  go = gsub('\\:', '-', go_df$go[str_to_upper(go_df$go_name) == str_to_upper(gene)])
  meta$expr = mat[go, meta$barcode]
  
  # extract data frame
  plot_df = data.frame()
  conditions = unique(meta$label)
  for (condition in conditions) {
    tmp_plot_df = meta %>% 
      filter(label == !!condition)
    fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.02, degree = 1)
    tmp_plot_df$expr = predict(fit, tmp_plot_df)
    plot_df %<>% rbind(tmp_plot_df)
  }
  title = go_df$go_name[gsub('\\:', '-', go_df$go) == go]
  
  range = range(plot_df$expr)
  brks = c(range[1] + 0.1 * diff(range),
           range[2] - 0.1 * diff(range))
  labels = c('min', 'max')
  
  pal = pals::ocean.thermal(100)
  pal = brewer.pal(9, 'RdGy') %>% rev
  pal = pals::ocean.solar(100)
  pal = nr_heat_blue_no_white %>% tail(-5)
  pal = nr_heat_blue_spatial
  p = plot_df %>%
    # mutate(label = ifelse(label == 'd1', 'Day 3', 'Day 7')) %>%
    # arrange(-expr) %>%
    ggplot(aes(x = x, y = y, fill = expr)) +
    rasterise(geom_point(size = 0.4,
                         shape = 21, stroke = 0, alpha = 1), dpi = 600) +
    scale_color_gradientn(colours = pal, name = 'Expression',
                          breaks = brks, labels = labels) +
    scale_fill_gradientn(colours = pal, name = 'Expression', 
                         breaks = brks, labels = labels) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    boxed_theme(size_lg = 6, size_sm = 5) +
    ggtitle(title) +
    coord_fixed() +
    theme(
      # aspect.ratio = 1,
      plot.title = element_text(size = 5, 
                                margin = margin(0,0,-5,0)),
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
      legend.justification = 'bottom',
      panel.background = element_rect(fill = 'grey96')
    ) +
    facet_wrap(~label)
  if (gene != last(genes_to_plot))
    p = p + theme(legend.position = 'none')
  p
  go_plots[[length(go_plots)+1]] = p
}
p2 = wrap_plots(go_plots, ncol = 2)
ggsave('fig/EFig9/GO-modules.pdf', p2, width = 10, height = 5.5,
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## d. lollipop plot: genes ####
###############################################################################-

# load results
dat0 = readRDS('data/real_data/DE_summaries/vespucci/Maniatis2019-seed=42-nsub=10-de=nebula_nbgmm.rds')
min_pval = min(dat0$p_val[dat0$p_val > 0])
dat0$p_val = vapply(dat0$p_val, function(x){max(min_pval, x)}, as.numeric(1))
dat0 %<>%
  mutate(up_reg = logFC > 0) %>%
  group_by(up_reg) %>%
  arrange(p_val, -abs(logFC)) %>%
  mutate(
    log_p_val = -log10(p_val)
  ) %>%
  slice(1:15)
dat0 %<>% 
  arrange(-logFC)
dat0$gene = factor(dat0$gene, levels=rev(as.character(dat0$gene)))

p3 = dat0 %>%
  ggplot(aes(x = gene, y = logFC)) +
  # facet_wrap(sign(logFC) ~ ., ncol = 1, scales = 'free') +
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
ggsave('fig/EFig9/lollipop-genes.pdf', p3, width = 8, height = 8, 
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## d. lollipop plot: GO modules ####
###############################################################################-

# load GO
go = get_ontology("data/GO/go.obo")

# load GLM results
dat0 = readRDS('data/real_data/DE_summaries/vespucci_GO/Maniatis2019-seed=42-nsub=10-de=glm.rds') %>% 
  filter(!is.na(p_val))
min_pval = min(dat0$p_val[dat0$p_val > 0])
dat0$p_val = vapply(dat0$p_val, function(x){max(min_pval, x)}, as.numeric(1))
dat0 %<>%
  mutate(up_reg = logFC > 0) %>%
  group_by(up_reg) %>%
  arrange(p_val, -abs(logFC)) %>%
  mutate(
    log_p_val = -log10(p_val)
  ) %>%
  slice(1:15)

# name GO terms
dat0 %<>%
  mutate(gene = chartr('-', ':', gene)) %>% 
  mutate(descr = go$name[gene]) 
dat0 %<>% 
  arrange(-logFC)
dat0$descr = factor(dat0$descr, levels=rev(as.character(dat0$descr)))

# plot
p4 = dat0 %>%
  ggplot(aes(x = descr, y = logFC)) +
  # facet_wrap(sign(logFC) ~ ., ncol = 1, scales = 'free') +
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
ggsave('fig/EFig9/lollipop-GO.pdf', p4, width = 18, height = 8, 
       units = 'cm', useDingbats = FALSE)
