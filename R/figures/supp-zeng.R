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
source('R/theme.R')

# Supp Fig 4a - registration
sc = readRDS('data/published_data/seurat/Zeng2023.rds')
rename_cell_types = c(
    'CA1' = 'CA1 excitatory neuron',
    'CA2' = 'CA2 excitatory neuron',
    'CA3' = 'CA3 excitatory neuron',
    'DG' = 'DG',
    'Astro' = 'Astrocyte',
    'Inh' = 'Inhibitory neuron',
    'Endo' = 'Endothelial cell',
    'CTX-Ex' = 'CTX excitatory neuron',
    'LHb' = 'LHb',
    'Micro' = 'Microglia',
    'Oligo' = 'Oligodendrocyte',
    'OPC' = 'OPC',
    'SMC' = 'SMC'
)

meta = sc@meta.data %>% 
    left_join(
        data.frame(
            'cell_type' = names(rename_cell_types),
            'cell_type_name' = unname(rename_cell_types)
        )
    ) %>%
    mutate(
       point_alpha = ifelse(marker == 'marker', 1, 1),
       # point_color = ifelse(marker == 'marker', cell_type, 'Background')
    )
# color_pal = c('grey50', nr_base_4) %>% setNames(c('Background', 'CA1', 'CA2', 'CA3', 'DG'))
meta$cell_type_name = factor(meta$cell_type_name, levels=unname(rename_cell_types))
color_pal = c(nr_base_11_light, "#FDB462", "#B3DE69") %>% setNames(levels(meta$cell_type_name))

p1 = meta %>% 
    ggplot(aes(x=ori_x, y=ori_y, fill=cell_type_name, color=cell_type_name, alpha=point_alpha)) +
    ggrastr::rasterise(
        geom_point(size = 0.2, shape = 21, stroke = 0), dpi = 300
    ) +
    coord_fixed() +
    boxed_theme() +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_blank(),
          axis.ticks.length.x = unit(0, 'lines'),
          axis.ticks.length.y = unit(0, 'lines'),
          legend.key.size = unit(0.4, 'lines'),
          legend.position = 'none',
          plot.title = element_text(size = 5)) +
    # facet_wrap(~replicate, nrow=1) +
    scale_alpha(guide = 'none')
# p1
p2 = meta %>% 
    ggplot(aes(x=x, y=y, fill=cell_type_name, color=cell_type_name, alpha=point_alpha)) +
    ggrastr::rasterise(
        geom_point(size = 0.2, shape = 21, stroke = 0), dpi = 300
    ) +
    coord_fixed() +
    boxed_theme() +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_blank(),
          axis.ticks.length.x = unit(0, 'lines'),
          axis.ticks.length.y = unit(0, 'lines'),
          legend.position = 'right',
          # legend.position = 'bottom',
          legend.justification = 'bottom',
          legend.title = element_blank(),
          legend.text = element_text(size=5),
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.18, 'lines'),
          legend.key.size = unit(0.4, 'lines'),
          plot.title = element_text(size = 5)) +
    # facet_wrap(~replicate, nrow=1) +
    guides(colour = guide_legend(override.aes = list(size=1))) +
    scale_alpha(guide = 'none')
# p2
p0 = wrap_plots(p1, p2, nrow=1)
# p0
ggsave('fig//EFig8/registration.pdf', p0, width=9, height = 6, units='cm')

# Supp Fig 4b - genes

ves_res = readRDS('data/published_data/vespucci/Zeng2023-seed=42-nsub=10.rds')
spatial_auc_res = ves_res$spatial_auc_result$aucs

dat0 = meta %>%
    left_join(spatial_auc_res) %>%
    filter(!is.na(auc))

fit = loess(auc ~ x * y, data = dat0, span = 0.02, degree = 1)
dat0$auc_fit = predict(fit, dat0)

range = range(dat0$auc_fit)
brks = c(range[1] + 0.1 * diff(range),
         range[2] - 0.1 * diff(range))
labels = format(range, digits = 2)

# new_color_pal = nr_heat_red_spatial
new_color_pal = auc_pal = pals::kovesi.linear_kryw_5_100_c67(100) %>% rev %>% tail(-5)
new_color_pal = cet_pal(100, name = 'l19') %>% rev()

fig7b = dat0 %>%
    arrange(auc_fit) %>%
    ggplot(aes(x = x, y = y, fill = auc_fit, color = auc_fit)) +
    ggrastr::rasterise(
        geom_point(size = 0.2, shape = 21, stroke = 0), dpi = 300
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
          legend.key.height = unit(0.18, 'lines'),
          plot.title = element_text(size = 5))
fig7b
ggsave('fig//EFig8/AUC.pdf', fig7b, width=6, height=6, units='cm')

meta = sc@meta.data %>%
    mutate(
        label = ifelse(grepl('control', label), 'Control', 'AD'),
        label = factor(label, levels=c('Control', 'AD'))
    )
expr = GetAssayData(sc, slot='counts')
colnames(expr) = gsub('-', '_', colnames(expr))

dat0 = meta %>% dplyr::select(barcode, x, y, label)
genes_to_plot = c(
    'C1QL3',
    'GNA14',
    'PROX1',
    'PLP1',
    'GRIN2B',
    'PTK2B',
    'HPCA',
    'NSF'
)

conditions = unique(dat0$label)
for (i in 1:length(genes_to_plot)) {
    if (i == 1) expr_plots = list()
    gene = genes_to_plot[i]
    dat0$expr = expr[str_to_title(gene), dat0$barcode]
    labels = c('min', 'max')
    
    plot_df = data.frame()
    for (condition in conditions){
        tmp_plot_df = dat0 %>% filter(label == !!condition)
        fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.02, degree = 1)
        tmp_plot_df$expr_fit = predict(fit, tmp_plot_df)
        plot_df %<>% rbind(tmp_plot_df)
    }
    
    range = range(plot_df$expr_fit)
    brks = c(range[1] + 0.1 * diff(range),
             range[2] - 0.1 * diff(range))
    # quan_range = quantile(plot_df$expr_fit, c(0, 0.9))
    
    color_pal = nr_heat_red_no_white %>% tail(-5)
    
    expr_plot = plot_df %>%
        mutate(expr_fit = winsorize(expr_fit, range)) %>%
        arrange(-expr_fit) %>%
        ggplot(aes(x = x, y = y, fill = expr_fit)) +
        rasterise(geom_point(size = 0.2,
                             shape = 21, stroke = 0, alpha = 1), dpi = 600) +
        scale_color_gradientn(colours = color_pal, name = 'Expression', 
                              # limits = quan_range,
                              breaks = brks,
                              # breaks = quan_range,
                              labels=labels) +
        scale_fill_gradientn(colours = color_pal, name = 'Expression', 
                             # limits = quan_range,
                             breaks = brks,
                             # breaks = quan_range,
                             labels=labels) +
        guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
        boxed_theme(size_lg = 5, size_sm = 5) +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        ggtitle(toupper(gene)) +
        theme(
            # aspect.ratio = 1,
            plot.title = element_text(size=5,margin = margin(0,0,-2,0)),
            plot.margin = unit(c(0.1, 0.05, 0.05, 0.5), "cm"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            # legend.position = 'none'
            legend.key.width = unit(0.18, 'lines'),
            legend.key.height = unit(0.18, 'lines'),
            legend.text = element_text(size = 5),
            legend.title = element_text(size = 5),
            legend.position = 'right',
            legend.justification = 'bottom'
        ) +
        facet_wrap(~label)
    if (i%%8 != 0){
        expr_plot = expr_plot + theme(legend.position = 'none')
    }
    expr_plots[[length(expr_plots)+1]] = expr_plot
}
fig7c = wrap_plots(expr_plots, nrow=4)
ggsave('fig//EFig8/genes.pdf', fig7c, width=10, height=11, units='cm')

source('R/theme.R')
dat0 = ves_res$de_feature_result
min_pval = min(dat0$p_val[dat0$p_val > 0])
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
ggsave('fig//EFig8/lollipop-genes.pdf', p3, width = 8, height = 8, 
       units = 'cm', useDingbats = FALSE)

meta = sc@meta.data %>%
    left_join(
        data.frame(
            'cell_type' = names(rename_cell_types),
            'cell_type_name' = unname(rename_cell_types)
        )
    ) %>%
    left_join(ves_res) %>%
    filter(!is.na(auc))

pairs = tidyr::crossing(celltype1 = unique(meta$cell_type_name),
                        celltype2 = unique(meta$cell_type_name)) %>% 
    filter(celltype1 != celltype2)

tests = pmap_dfr(pairs, function(...) {
    current = tibble(...)
    print(current)
    
    vec1 = filter(meta, cell_type_name == current$celltype1) %>% 
        arrange(barcode) %>% pull(auc)
    vec2 = filter(meta, cell_type_name == current$celltype2) %>% 
        arrange(barcode) %>% pull(auc)
    t = t.test(vec1, vec2, alternative = 'g')$p.value
    w = wilcox.test(vec1, vec2, alternative = 'g')$p.value
    data.frame(test = c('t-test', 'wilcox'),
               pval = c(t, w)) %>% 
        cbind(current, .)
})

# plot median heatmap
delta = pmap_dfr(pairs, function(...) {
    current = tibble(...)
    vec1 = filter(meta, cell_type_name == current$celltype1) %>% 
        arrange(barcode) %>% pull(auc)
    vec2 = filter(meta, cell_type_name == current$celltype2) %>% 
        arrange(barcode) %>% pull(auc)
    median = median(vec1 - vec2)
    mutate(current, delta = median)
})
lvls = with(meta, reorder(cell_type_name, auc, stats::median)) %>% levels()
range = range(delta$delta)

labels = tests %>% 
    filter(test == 'wilcox') %>% 
    mutate(celltype1 = factor(celltype1, levels = lvls),
           celltype2 = factor(celltype2, levels = lvls),
           lab = ifelse(pval < 0.001/2, '***',
                        ifelse(pval < 0.01/2, '**',
                               ifelse(pval < 0.05/2, '*', ''))))
p1 = delta %>% 
    mutate(celltype1 = factor(celltype1, levels = lvls),
           celltype2 = factor(celltype2, levels = lvls)) %>% 
    ggplot(aes(x = celltype2, y = celltype1)) +
    geom_tile(color = 'white', aes(fill = delta)) +
    geom_text(data = labels, size = 1, aes(label = lab)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                           name = expression(Delta~AUC),
                           breaks = range,
                           labels = format(range, digits = 2)) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black', title.position='top')) +
    coord_fixed() +
    boxed_theme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.key.height = unit(0.18, 'lines'),
          legend.key.width = unit(0.25, 'lines'),
          legend.position = 'bottom',
          legend.justification = 'right')
ggsave(paste0('fig//EFig8/celltype-auc-delta-wilcox.pdf'), p1, width=6, height=6, units='cm')

