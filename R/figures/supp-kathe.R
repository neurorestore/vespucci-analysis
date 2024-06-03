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

# Supp Fig 6a marker genes 
sc = readRDS('data/published_data/seurat/Kathe2022.rds')
meta = sc@meta.data %>%
    mutate(
        delta_x = ori_x - x,
        x = ori_x
    )

lamina_dat = readRDS('data/published_data/raw_data/Kathe2022/medres.rds') %>%
    dplyr::rename(barcode = barcodes) %>%
    mutate(
        barcode = gsub('-', '_', barcode)
    )
all(colnames(meta$barcode) %in% lamina_dat$barcode)

meta %<>% 
    left_join(
        lamina_dat
    )

meta$replicate = vapply(meta$replicate, function(x){
    case_when(
        x == 'M2' ~ 'EES-REHAB Rep 1',
        x == 'M3' ~ 'SCI Rep 1',
        x == 'M5' ~ 'SCI Rep 2',
        x == 'M7' ~ 'EES REHAB Rep 2',
        x == 'M9' ~ 'EES-REHAB Rep 3',
        x == 'M10' ~ 'SCI Rep 3'
    )
}, as.character(1))

meta$medres_clean = factor(vapply(as.character(meta$medres), function(x){
    case_when(
        x == 'lamina-1-4' ~ 'Lamina layers 1-4',
        x == 'lamina-5' ~ 'Lamina layer 5',
        x == 'intermediate-dorsal' ~ 'Intermediate dorsal',
        x == 'intermediate-ventral' ~ 'Intermediate ventral',
        x == 'ventral' ~ 'Ventral'
    )
}, as.character(1)), levels=c(
    'Lamina layers 1-4','Lamina layer 5','Intermediate dorsal','Intermediate ventral','Ventral'
))

color_pal = nr_base_6[2:6] %>% setNames(levels(meta$medres_clean[!is.na(meta$medres_clean)]))

fig6a = meta %>% 
    ggplot(aes(x=x, y=y, fill=medres_clean, color=medres_clean)) +
    ggrastr::rasterise(
        geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
    ) +
    coord_fixed() +
    boxed_theme() +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y = element_blank(),
          # strip.text.x = element_blank(),
          axis.ticks.length.x = unit(0, 'lines'),
          axis.ticks.length.y = unit(0, 'lines'),
          # legend.position = 'none',
          legend.position = 'right',
          legend.justification = 'bottom',
          legend.title = element_blank(),
          legend.key.size = unit(0.4, 'lines'),
          plot.title = element_text(size = 5)
          ) 
    # guides(colour = guide_legend(override.aes = list(size=2)))
# fig6a
ggsave('fig/final/Efig7/lamina.pdf', fig6a, width=7, height = 8, units='cm')

ves_res = readRDS('data/published_data/vespucci/Kathe2022-seed=42-nsub=10.rds')
spatial_auc_res = ves_res$spatial_auc_result$aucs

comparison = 'EES_REHAB-SCI'

dat0 = meta %>%
    left_join(spatial_auc_res) %>%
    filter(!is.na(auc))

fit = loess(auc ~ x * y, data = dat0, span = 0.02, degree = 1)
dat0$auc_fit = predict(fit, dat0)

range = range(dat0$auc_fit)
brks = c(range[1] + 0.1 * diff(range),
         range[2] - 0.1 * diff(range))
labels = format(range, digits = 2)

new_color_pal = cet_pal(100, name = 'l19') %>% rev()

fig6b = dat0 %>%
    arrange(auc_fit) %>%
    ggplot(aes(x = x, y = y, fill = auc_fit, color = auc_fit)) +
    ggrastr::rasterise(
        geom_point(size = 0.7, shape = 21, stroke = 0), dpi = 300
    ) +
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
# p1
ggsave('fig/final/Efig7/AUC.pdf', fig6b, width=6, height=7, units='cm')

expr = GetAssayData(sc, slot='counts')
colnames(expr) = gsub('-', '_', colnames(expr))

expr %<>% NormalizeData()

dat0 = meta %>% dplyr::select(barcode, x, y, label)
genes_to_plot = c(
    'CAMK2N1',
    'FGF1',
    'SCG2',
    'AHI1',
    'ZCCHC12',
    'SYT1'
)

sum(!str_to_title(genes_to_plot) %in% rownames(expr))

conditions = unique(dat0$label)
for (i in 1:length(genes_to_plot)) {
    if (i == 1) expr_plots = list()
    gene = genes_to_plot[i]
    dat0$expr = expr[str_to_title(gene), dat0$barcode]
    expr_range = c()
    for (condition in conditions){
        tmp_plot_df = dat0 %>% dplyr::filter(label == !!condition)
        fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.02, degree = 1)
        tmp_plot_df$expr = predict(fit, tmp_plot_df)
        expr_range = c(expr_range, tmp_plot_df$expr)
    }
    
    full_range = range(expr_range)
    brks = range(full_range)
    # labels = format(full_range, digits = 3)
    labels = c('min', 'max')
    
    plot_df = data.frame()
    for (condition in conditions){
        tmp_plot_df = dat0 %>% filter(label == !!condition)
        fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.02, degree = 1)
        tmp_plot_df$expr_fit = predict(fit, tmp_plot_df)
        tmp_plot_df$expr_fit = winsorize(tmp_plot_df$expr_fit, full_range)
        plot_df %<>% rbind(tmp_plot_df)
    }
    # color_pal = colorRampPalette(c("lightyellow", "yellow", "purple"))(10)
    color_pal = nr_heat_red_no_white %>% tail(-5)
    
    plot_df$label = gsub('_', ' ', plot_df$label)
    plot_df$label = factor(plot_df$label, levels=c('SCI', 'EES REHAB'), 
                           labels = c(expression(SCI), expression('EES'^REHAB)))
    plot_df$alpha = 1
    plot_df$alpha[plot_df$x < 250] = 1
    
    expr_plot = plot_df %>%
        arrange(-expr_fit) %>%
        ggplot(aes(x = x, y = y, fill = expr_fit, alpha=alpha)) +
        rasterise(geom_point(size = 0.5,
                             shape = 21, stroke = 0), dpi = 600) +
        scale_color_gradientn(colours = color_pal, name = 'Expression', breaks = brks, labels=labels) +
        scale_fill_gradientn(colours = color_pal, name = 'Expression', breaks = brks, labels=labels) +
        guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
        scale_alpha(guide = 'none') +
        boxed_theme(size_lg = 6, size_sm = 5) +
        ggtitle(toupper(gene)) +
        coord_fixed() +
        theme(
            # aspect.ratio = 1,
            plot.title = element_text(size=5,margin = margin(0,0,-5,0)),
            plot.margin = unit(c(0.05, 0.05, 0.05, 0.5), "cm"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            # legend.position = 'none'
            legend.key.width = unit(0.18, 'lines'),
            legend.key.height = unit(0.14, 'lines'),
            legend.text = element_text(size = 3.5),
            legend.title = element_text(size = 3.5),
            legend.position = 'right',
            legend.justification = 'bottom'
        ) +
        # geom_vline(xintercept = 250, linetype = 'dashed', size=0.15, color='black') +
        facet_wrap(~label, labeller = label_parsed)
    if (i%%6 != 0){
        expr_plot = expr_plot + theme(legend.position = 'none')
    }
    expr_plots[[length(expr_plots)+1]] = expr_plot
}
fig6c = wrap_plots(expr_plots, nrow=2)
ggsave(paste0('fig/final/Efig7/genes.pdf'), fig6c, width=18, height=6, units='cm')

## # Supp Fig 4c - GO
sc = readRDS('data/published_data/seurat_GO/Kathe2022.rds')
meta = sc@meta.data %>%
    mutate(
        x = ori_x
    )

expr = GetAssayData(sc, slot='counts')
colnames(expr) = gsub('-', '_', colnames(expr))

dat0 = meta %>% dplyr::select(barcode, x, y, label)
go_df = readRDS('data/metadata/go_names.rds')

genes_to_plot = c(
    'NEUROPEPTIDE HORMONE ACTIVITY',
    'GLUTAMATE DECARBOXYLASE ACTIVITY',
    'SYNCHRONOUS NEUROTRANSMITTER SECRETION',
    'RESPIRATORY CHAIN COMPLEX IV',
    'SENSORY PERCEPTION',
    'FAST, CALCIUM ION-DEPENDENT EXOCYTOSIS OF NEUROTRANSMITTER'
)
conditions = unique(dat0$label)
for (i in 1:length(genes_to_plot)) {
    if (i == 1) expr_plots = list()
    gene = genes_to_plot[i]
    
    go = gsub('\\:', '-', go_df$go[str_to_upper(go_df$go_name) == str_to_upper(gene)])
    dat0$expr = expr[go, dat0$barcode]
    expr_range = c()
    for (condition in conditions){
        tmp_plot_df = dat0 %>% dplyr::filter(label == !!condition)
        fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.02, degree = 1)
        tmp_plot_df$expr = predict(fit, tmp_plot_df)
        expr_range = c(expr_range, tmp_plot_df$expr)
    }
    
    full_range = range(expr_range)
    if (gene == 'RESPIRATORY CHAIN COMPLEX IV') {
        full_range = quantile(tmp_plot_df$expr, c(0.1, 0.99))   
    }
    brks = quantile(full_range)[-1] %>% unname()
    # labels = format(full_range, digits = 3)
    labels = c('min', '', '', 'max')
    
    plot_df = data.frame()
    for (condition in conditions){
        tmp_plot_df = dat0 %>% filter(label == !!condition)
        fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.02, degree = 1)
        tmp_plot_df$expr_fit = predict(fit, tmp_plot_df)
        tmp_plot_df$expr_fit = winsorize(tmp_plot_df$expr_fit, full_range)
        plot_df %<>% rbind(tmp_plot_df)
    }
    # color_pal = colorRampPalette(c("yellow", "purple"))(10)
    color_pal = nr_heat_blue_spatial
    go_name_clean = go_df$go_name[str_to_upper(go_df$go_name) == str_to_upper(gene)]
    
    plot_df$label = gsub('_', ' ', plot_df$label)
    plot_df$label = factor(plot_df$label, levels=c('SCI', 'EES REHAB'), 
                           labels = c(expression(SCI), expression('EES'^REHAB)))
    
    expr_plot = plot_df %>%
        arrange(-expr_fit) %>%
        ggplot(aes(x = x, y = y, fill = expr_fit)) +
        rasterise(geom_point(size = 0.5,
                             shape = 21, stroke = 0), dpi = 600) +
        scale_color_gradientn(colours = color_pal, name = 'GO module\nscore', breaks = brks, labels=labels) +
        scale_fill_gradientn(colours = color_pal, name = 'GO module\nscore', breaks = brks, labels=labels) +
        guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
        scale_alpha(guide = 'none') +
        boxed_theme(size_lg = 6, size_sm = 5) +
        ggtitle(go_name_clean) +
        coord_fixed() +
        theme(
            # aspect.ratio = 1,
            plot.title = element_text(size=5,margin = margin(0,0,-5,0)),
            plot.margin = unit(c(0.05, 0.05, 0.05, 0.5), "cm"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            # legend.position = 'none'
            legend.key.width = unit(0.18, 'lines'),
            legend.key.height = unit(0.18, 'lines'),
            legend.text = element_text(size = 3.5),
            legend.title = element_text(size = 3.5),
            legend.position = 'right',
            legend.justification = 'bottom'
        ) +
        # geom_vline(xintercept = 250, linetype = 'dashed', size=0.15, color='black') +
        facet_wrap(~label, labeller = label_parsed)
    if (i%%6 != 0){
        expr_plot = expr_plot + theme(legend.position = 'none')
    }
    expr_plots[[length(expr_plots)+1]] = expr_plot
}
fig6d = wrap_plots(expr_plots, nrow=2)
ggsave('fig/final/Efig7/GO.pdf', fig6d, width=18, height=6, units='cm')

# Supp Fig 6e
source('R/theme.R')
dat0 = ves_res$$de_feature_result
min_pval = min(dat0$p_val[dat0$p_val > 0])
dat0$p_val = vapply(dat0$p_val, function(x){max(min_pval, x)}, as.numeric(1))
dat0 %<>%
    mutate(up_reg = effect_size > 0) %>%
    group_by(up_reg) %>%
    arrange(p_val, -abs effect_size)) %>%
    mutate(
        log_p_val = -log10(p_val)
    ) %>%
    slice(1:15)

dat0 %<>% 
    arrange( effect_size)
dat0$gene = factor(dat0$gene, levels=rev(as.character(dat0$gene)))

fig6e = dat0 %>%
    ggplot(aes(x = gene, y = effect_size)) +
    # facet_wrap(sign effect_size) ~ ., ncol = 1, scales = 'free') +
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
fig6e
ggsave('fig/final/Efig7/lollipop-genes.pdf', fig6e, width=8, height=8, units='cm')


# Supp Fig 6f
source('R/theme.R')
dat0 = readRDS('data/published_data/vespucci_GO/Kathe2022-seed=42-nsub=10.rds')$de_feature_result
dat0 %<>%
    filter(!is.na(p_val))
min_pval = min(dat0$p_val[dat0$p_val > 0])
dat0$p_val = vapply(dat0$p_val, function(x){max(min_pval, x)}, as.numeric(1))
dat0 %<>%
    mutate(up_reg = effect_size > 0) %>%
    group_by(up_reg) %>%
    arrange(p_val, -abs effect_size)) %>%
    mutate(
        log_p_val = -log10(p_val)
    ) %>%
    slice(1:15)

go_df = readRDS('data/metadata/go_names.rds') %>%
    mutate(gene = gsub('\\:', '-', go))
dat0 %<>% left_join(
    go_df, by='gene'
) %>%
    arrange(p_val, -abs effect_size))

dat0 %<>% 
    arrange( effect_size)
dat0$go_name = factor(dat0$go_name, levels=rev(as.character(dat0$go_name)))

fig6f= dat0 %>%
    ggplot(aes(x = go_name, y = effect_size, color=log_p_val, fill=log_p_val)) +
    geom_hline(aes(yintercept = 0), linetype = "dotted", size = 0.3, color='grey80') +
    geom_segment(aes(xend = go_name, yend = 0), color='grey80', alpha=0.5) +
    geom_point(shape = 21, stroke = 0.2, size = 0.8, color = 'black', fill='grey80') +
    scale_y_continuous(expression(beta)) +
    coord_flip() +
    boxed_theme() +
    scale_color_gradientn(
        colors = color_pal,
        name = expression("-log"[10]~"(P)"),
        breaks = brks,
        labels = labels
    ) +
    scale_fill_gradientn(
        colors = color_pal,
        name = expression("-log"[10]~"(P)"),
        breaks = brks,
        labels = labels
    ) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    theme(
        axis.title.y = element_blank(),
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.2, 'lines'),
        legend.key.height = unit(0.2, 'lines'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
    )
fig6f
ggsave('fig/final/Efig7/lollipop-GO.pdf', fig6f, width=10, height=8, units='cm')
