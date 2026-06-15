# setwd('C:/Users/teo/Documents/EPFL/projects/vespucci/')
setwd('~/git/vespucci-analysis/')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(pROC)
library(PRROC)
library(lawstat)
library(nparcomp)
library(gridExtra)
library(cetcolor)
source('R/theme.R')

#############################################################################-
# Calcagno registration
#############################################################################-
sc = readRDS('data/real_data/seurat/Calcagno2022.rds')
meta = sc@meta.data %>%
    mutate(
        barcode = gsub('-', '_', barcode),
        label = timepoint
    )

meta %>% dplyr::select(annotation, annotation_name) %>% table()
meta$annotation_name = factor(meta$annotation_name, levels = c('Remote zone', 'Border zone 1', 'Border zone 2', 'Infarction zone'))

color_pal = c('#746FB4', '#AB4567', '#E3201D', '#BEBEBE') %>% set_names(levels(meta$annotation_name))
p6_1 = meta %>%
    ggplot(aes(x=ori_x,y=ori_y, fill=annotation_name, color=annotation_name)) +
    ggrastr::rasterise(
        geom_point(size = 0.3, shape = 21, stroke = 0), dpi = 600
    ) +
    boxed_theme() +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    theme(
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.5), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0, 'lines'),
        axis.ticks.length.y = unit(0, 'lines'),
        legend.position = 'none',
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) + 
    facet_wrap(~ replicate, nrow = 1, scales='free')
p6_2 = meta %>%
    ggplot(aes(x=x,y=y, fill=annotation_name, color=annotation_name)) +
    ggrastr::rasterise(
        geom_point(size = 0.3, shape = 21, stroke = 0), dpi = 600
    ) +
    boxed_theme() +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    theme(
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.5), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0, 'lines'),
        axis.ticks.length.y = unit(0, 'lines'),
        legend.position = 'bottom',
        legend.justification = 'bottom',
        legend.text = element_text(size=5),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.spacing.x = unit(0.01, 'cm')
    ) + 
    guides(colour = guide_legend(override.aes = list(size=1))) +
    facet_wrap(~ replicate, nrow = 1, scales='free')
p6 = wrap_plots(p6_1, p6_2, nrow=2)
ggsave('fig/final/Fig1/calcagno_registration.pdf', p6, width = 7, height = 3.2, units='cm')


#############################################################################-
# Calcagno auc
#############################################################################-
auc_res = readRDS('data/real_data/vespucci/Calcagno2022-seed=42-nsub=10.rds')[[1]]$aucs
meta %<>% left_join(auc_res)

# interpolate in 2D
fit = loess(auc ~ x * y, data = meta, span = 0.015)
meta$auc_fit = predict(fit, meta)

range = range(meta$auc_fit)
brks = c(range[1] + 0.1 * diff(range),
         range[2] - 0.1 * diff(range))
labels = format(range, digits = 2)
labels = c(paste0(labels[1], ' '),
           paste0(' ', labels[2]))

auc_pal = cet_pal(100, name = 'l19') %>% rev()

p7 = meta %>%
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
    guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE),
           color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
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
          legend.justification = 'bottom',
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.18, 'lines'),
          # plot.title = element_text(size = 5)
          plot.title = element_blank()
    )
p7
ggsave(paste0("fig/final/Fig1/calcagno_AUC.pdf"), p7, width = 2, height = 2.7, units = "cm", useDingbats = FALSE)


###############################################################################-
## Calcagno genes ####
###############################################################################-
sc = readRDS('data/real_data/seurat/Calcagno2022.rds')
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
    'SPP1', 'SFRP2'
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
            plot.title = element_text(size = 5, margin = margin(0,0,-2,0)),
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
            legend.position = 'bottom',
            legend.justification = 'right'
        ) +
        facet_wrap(~label) +
        coord_fixed()
    if (idx != 2)
        p = p + theme(legend.position = 'none')
    p
    expr_plots[[idx]] = p
}
p8 = wrap_plots(expr_plots, nrow = 1)
ggsave('fig/final/Fig1/calcagno_genes.pdf', p8, width = 5, height = 2.7,
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## Calcagno GO terms ####
###############################################################################-

# load GO
sc = readRDS('data/real_data/seurat_GO/DE/Calcagno2022.rds')[[1]]
go_df = readRDS('data/metadata/go_names.rds')
meta = sc@meta.data %>%
    mutate(barcode = gsub('-', '_', barcode))
mat = GetAssayData(sc, slot='counts')
gos_to_plot = c(
    'response to injury involved in regulation of muscle adaptation',
    'negative regulation of antigen processing and presentation of peptide or polysaccharide antigen via MHC class II'
)
gos = go_df %>% filter(go_name %in% gos_to_plot) %>%
    pull(go) %>%
    gsub('\\:', '-', .)
mat = mat[gos,]

# iterate through GO terms
go_plots = list()
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
        scale_color_gradientn(colours = pal, name = 'GO module\nscore',
                              breaks = brks, labels = labels) +
        scale_fill_gradientn(colours = pal, name = 'GO module\nscore', 
                             breaks = brks, labels = labels) +
        guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
        boxed_theme(size_lg = 6, size_sm = 5) +
        ggtitle(title) +
        theme(
            aspect.ratio = 1,
            plot.title = element_text(size = 5, margin = margin(0,0,-2,0)),
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
            legend.position = 'bottom',
            legend.justification = 'right'
        ) +
        facet_wrap(~label) + 
        coord_fixed()
    if (idx != 2)
        p = p + theme(legend.position = 'none')
    p
    go_plots[[idx]] = p
}
p9 = wrap_plots(go_plots, nrow = 1)
ggsave('fig/final/Fig1/calcagno_GO-modules.pdf', p9, width = 5, height = 3.5,
       units = 'cm', useDingbats = FALSE)


