setwd('~/git/vespucci-analysis/')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(lawstat)
library(nparcomp)
library(ontologyIndex)
library(cetcolor)
source('R/theme.R')


###############################################################################-
## a. RCTD ####
###############################################################################-

meta = readRDS('data/real_data/seurat/Calcagno2022.rds')@meta.data %>%
    mutate(
        barcode = gsub('-', '_', barcode),
        label = timepoint
    )
rctd_res = readRDS('data/real_data/rctd/Calcagno2022/rctd-summary-res.rds') %>%
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
ggsave('fig/EFig8/rctd-registration.pdf', p0, width = 11, height = 4, units='cm')

###############################################################################-
## b. genes ####
###############################################################################-

# read Calcagno dataset
sc = readRDS('data/real_data/seurat/Calcagno2022.rds')
meta = sc@meta.data %>%
  mutate(
    barcode = gsub('-', '_', barcode),
    label = timepoint
  )
expr = GetAssayData(sc, slot = 'counts')
colnames(expr) = gsub('-', '_', colnames(expr))

# normalize
normalize_data = TRUE
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
ggsave('fig/EFig8/genes.pdf', p1, width = 10, height = 4.5,
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## c. GO terms ####
###############################################################################-

# load GO
# go = get_ontology("data/GO/go.obo")
# 
# # read GO sub-matrix
# mat = readRDS('data/real_data/seurat_GO/Calcagno2022_supp_subset_GO.rds')
# meta = data.frame(barcode = colnames(mat)) %>% 
#   separate(barcode, into = c('condition', 'label', 'replicate', 'x'),
#            sep = '_', remove = FALSE) %>% 
#   dplyr::select(-x)
# # merge in coordinates from single-cell object
# coords = sc@meta.data %>% 
#   dplyr::select(barcode, x, y) %>% 
#   mutate(barcode = chartr('-', '_', barcode))
# meta %<>% left_join(coords, by = 'barcode')

sc = readRDS('data/real_data/seurat_GO/DE/Calcagno2022.rds')[[1]]
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
# go_plots = list()
go_plots = expr_plots
for (idx in seq_len(nrow(mat))) {
  # go_term = rownames(mat)[idx] %>% 
  #   chartr('-', ':', .)
  # descr = go$name[go_term]
  # title = paste0(go_term, '\n', descr)
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
# ggsave('fig/EFig8/GO-modules.pdf', p2, width = 12, height = 3,
       # units = 'cm', useDingbats = FALSE)
ggsave('fig/EFig8/genes-GO-modules.pdf', p2, width = 16, height = 5,
       units = 'cm', useDingbats = FALSE)


###############################################################################-
## d. lollipop plot: genes ####
###############################################################################-

# load results
dat0 = readRDS('data/real_data/DE_summaries/vespucci/Calcagno2022-seed=42-nsub=10-de=nebula_nbgmm.rds')
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
ggsave('fig/EFig8/lollipop-genes.pdf', p3, width = 8, height = 8, 
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## e. lollipop plot: GO modules ####
###############################################################################-

# load GO
go = get_ontology("data/GO/go.obo")

# load GLM results
dat0 = readRDS('data/real_data/DE_summaries/vespucci_GO/Calcagno2022-seed=42-nsub=10-de=glm.rds') %>% 
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
ggsave('fig/EFig8/lollipop-GO.pdf', p4, width = 18, height = 8, 
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## unregister AUCs
###############################################################################-

# load GO
# meta = readRDS('data/real_data/meta/Calcagno2022.rds') %>% mutate(barcode = gsub('-','_',barcode))
# aucs = readRDS('data/real_data/vespucci/Calcagno2022-seed=42-nsub=10.rds')[[1]]$aucs 
# meta %<>% left_join(aucs)
sc = readRDS('data/real_data/seurat/Calcagno2022.rds')
meta = sc@meta.data %>%
  mutate(
    barcode = gsub('-', '_', barcode),
    label = timepoint
  )
expr = GetAssayData(sc, slot = 'counts')
colnames(expr) = gsub('-', '_', colnames(expr))
expr %<>% NormalizeData()
aucs = readRDS('data/real_data/vespucci/Calcagno2022-seed=42-nsub=10.rds')[[1]]$spatial_auc_result$aucs 
meta %<>% left_join(aucs)
# p1_3

new_meta = map_df(unique(meta$replicate), function(rep){
  tmp_meta = meta %>% filter(replicate == rep)
  fit = loess(auc ~ ori_x * ori_y, data = tmp_meta, span = 0.015)
  tmp_meta$auc_fit = predict(fit, tmp_meta)
  tmp_meta
})

range = range(new_meta$auc_fit)
brks = c(range[1] + 0.1 * diff(range), range[2] - 0.1 * diff(range))
labels = format(range, digits = 2)
labels = c(paste0(labels[1], ' '), paste0(' ', labels[2]))
auc_pal = cet_pal(100, name = 'l19') %>% rev()
p5_1 = new_meta %>%
    mutate(
      replicate_clean = paste0(ifelse(label == 'd1', 'Day 3 ', 'Day 7 '), str_to_title(gsub('.*_', '', replicate)))
    ) %>%
    # arrange(auc_fit) %>%
    ggplot(aes(x = ori_x, y = ori_y, fill = auc_fit, color = auc_fit)) +
    ggtitle('Vespucci') +
    ggrastr::rasterise(
        geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 600
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = auc_pal, name = 'AUC   ', labels = labels, limits = range, breaks = brks) +
    scale_color_gradientn(colours = auc_pal, name = 'AUC   ', labels = labels, limits = range, breaks = brks) +
    guides(fill = guide_colorbar(ticks.colour = NA, frame.colour = NA, size=0.1, title.position='top')) +
    # coord_fixed() +
    boxed_theme(size_lg = 5, size_sm = 5) +
    theme(
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
		legend.position = 'bottom',
		legend.justification = 'bottom',
		legend.key.width = unit(0.18, 'lines'),
		legend.key.height = unit(0.18, 'lines'),
        # plot.title = element_blank()
    ) +
    facet_wrap(~ replicate_clean, nrow = 2, scales='free')
# p5
plot_list = list(p5_1)
for (gene in c('Spp1', 'Sfrp2')) {
	new_meta$expr = expr[gene, new_meta$barcode]
	replicates = unique(new_meta$replicate)
	plot_df = data.frame()
	for (replicate in replicates) {
		tmp_plot_df = new_meta %>% dplyr::filter(replicate == !!replicate)
		fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.02, degree = 1)
		tmp_plot_df$expr = predict(fit, tmp_plot_df)
		plot_df %<>% rbind(tmp_plot_df)
	}
	
	range = range(plot_df$expr)
	brks = c(range[1] + 0.1 * diff(range), range[2] - 0.1 * diff(range))
	labels = c('min', 'max')
	pal = nr_heat_red_no_white %>% tail(-5)
	p5_2 = plot_df %>%
		mutate(replicate_clean = paste0(ifelse(label == 'd1', 'Day 3 ', 'Day 7 '), str_to_title(gsub('.*_', '', replicate)))) %>%
		# arrange(auc_fit) %>%
		ggplot(aes(x = ori_x, y = ori_y, fill = expr, color = expr)) +
		ggtitle('Vespucci') +
		ggrastr::rasterise(
			geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 600
		) +
		scale_y_continuous(expand = c(0, 0)) +
		scale_x_continuous(expand = c(0, 0)) +
		scale_fill_gradientn(colours = pal, name = 'Expr   ', labels = labels, limits = range, breaks = brks) +
		scale_color_gradientn(colours = pal, name = 'Expr   ', labels = labels, limits = range, breaks = brks) +
		guides(fill = guide_colorbar(ticks.colour = NA, frame.colour = NA, size=0.1, title.position='top')) +
		# coord_fixed() +
		boxed_theme(size_lg = 5, size_sm = 5) +
		theme(
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
			legend.position = 'bottom',
			legend.justification = 'bottom',
			legend.key.width = unit(0.18, 'lines'),
			legend.key.height = unit(0.18, 'lines'),
			# plot.title = element_blank()
		) +
		facet_wrap(~ replicate_clean, nrow = 2, scales='free')
	plot_list[[length(plot_list)+1]] = p5_2
}
p5 = wrap_plots(plot_list, nrow=1)
ggsave('fig/EFig8/unregister-aucs.pdf', p5, width=13.5, height=6, units='cm')

meta = readRDS('data/real_data/meta/Calcagno2022.rds') %>% 
    mutate(
        replicate_clean = paste0(ifelse(label == 'd1', 'Day 3 ', 'Day 7 '), str_to_title(gsub('.*_', '', replicate)))
    )
replicates = unique(meta$replicate_clean)
ves_files = list.files('data/real_data/vespucci_leave_one_out/', pattern='remove', full.names=T)

plot_list = list()
for (ves_file in ves_files) {
    ves_res = readRDS(ves_file)
    meta0 = meta %>% 
        mutate(barcode = gsub('-','_',barcode)) %>% 
        left_join(ves_res[[1]]$spatial_auc_result$aucs) %>%
        filter(!is.na(auc))
    replicate_removed = replicates[!replicates %in% meta0$replicate_clean]
    stopifnot(length(replicate_removed) == 1)
    
    # interpolate in 2D
    fit = loess(auc ~ x * y, data = meta0, span = 0.015)
    meta0$auc_fit = predict(fit, meta0)
    
    range = range(meta0$auc_fit)
    brks = c(range[1] + 0.1 * diff(range),
             range[2] - 0.1 * diff(range))
    labels = format(range, digits = 2)
    labels = c(paste0(labels[1], ' '),
               paste0(' ', labels[2]))
    auc_pal = cet_pal(100, name = 'l19') %>% rev()
    p1 = meta0 %>%
        # arrange(auc_fit) %>%
        ggplot(aes(x = x, y = y, fill = auc_fit, color = auc_fit)) +
        ggtitle(paste0(replicate_removed, '\nremoved')) +
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
            legend.position = 'bottom',
            legend.justification = 'bottom',
            legend.key.width = unit(0.18, 'lines'),
            legend.key.height = unit(0.18, 'lines'),
            # plot.title = element_blank()
        )
    # p1
    plot_list[[length(plot_list)+1]] = p1
}
p6 = wrap_plots(plot_list, nrow=2)
p6
ggsave('fig/EFig8/leave-rep-out-aucs.pdf', p6, width=5, height=7, units='cm')

pairs = tidyr::crossing(
    ves_file1 = ves_files,
    ves_file2 = ves_files
) %>% filter(
    ves_file1 != ves_file2
) %>% mutate(
    ves1_remove = str_to_title(gsub('d7_', 'Day 7 ', gsub('d1_', 'Day 3 ', gsub('MI_', '', gsub('-seed=42.rds', '', gsub('.*remove=', '', ves_file1)))))),
    ves2_remove = str_to_title(gsub('d7_', 'Day 7 ', gsub('d1_', 'Day 3 ', gsub('MI_', '', gsub('-seed=42.rds', '', gsub('.*remove=', '', ves_file2))))))
)

cor_df = map_df(1:nrow(pairs), function(i){
    tmp_row = pairs[i,]
    de_res = readRDS(tmp_row$ves_file1)[[1]]$de_feature_result %>% dplyr::select(feature, p_val) %>% dplyr::rename(ves1_pval = p_val) %>% inner_join(readRDS(tmp_row$ves_file2)[[1]]$de_feature_result %>% dplyr::select(feature, p_val) %>% dplyr::rename(ves2_pval = p_val))
    tmp_row %>% dplyr::select(ves1_remove, ves2_remove) %>% mutate(pearson_cor = cor(de_res$ves1_pval, de_res$ves2_pval, method='pearson', use='complete.obs'), spearman_cor = cor(de_res$ves1_pval, de_res$ves2_pval, method='spearman', use='complete.obs'))
})

cor_df$val = cor_df$spearman_cor
range = range(cor_df$val)
brks = c(range[1] + 0.1 * diff(range),
         range[2] - 0.1 * diff(range))

p_out = cor_df %>% 
    ggplot(aes(x = ves1_remove, y = ves2_remove)) +
    geom_tile(color = 'white', aes(fill = val)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",name = 'Corr.',breaks = brks, labels = format(range, digits = 2)) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black', title.position = 'top')) +
    coord_fixed() +
    boxed_theme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.15, 'lines'),
          legend.position = 'bottom',
          legend.justification = 'right')
p_out

ggsave('fig/EFig8/leave-one-out-cor.pdf', p_out, width=4, height=4, units='cm')
