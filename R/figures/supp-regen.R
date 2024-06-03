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
library(ungeviz)
source('R/theme.R')

source("R/kinematics/kinematics_parsing.R")
dat = readRDS('data/kinematics/final_data.rds')
meta_cols = colnames(dat)[1:4]

p1_col = 'drag_duration_percentage'
p2_col = 'step_height_foot'
p3_col = 'ampl_elev_tight'

p1_name = "% Drag"
p2_name = "Step height"
p3_name = "Thigh oscillation"

p1_lab = "% Drag (% total time)"
p2_lab = "Step height of foot (cm)"
p3_lab = "Thigh oscillation (degrees)"

factor_colors = c(
    'young' = 'grey70',
    'old' = 'grey10',
    'treated' = nr_tree_red[70]
)

kin = parse_kinematics(dat,
                       # Meta data info
                       meta_cols = meta_cols,
                       label_col = 'condition',
                       replicate_col = 'replicate',
                       # Variables to plot
                       p1_col = p1_col,
                       p2_col = p2_col,
                       p3_col = p3_col,
                       pc_flip = F,
                       plot_data = T,
                       factor_colors = factor_colors)

p0 = plot_kinematics(kin, 
                     # Names of these variables for the loadings plot
                     p1_name = p1_name,
                     p2_name = p2_name,
                     p3_name = p3_name,
                     # Axis labels
                     p1_lab = p1_lab,
                     p2_lab = p2_lab,
                     p3_lab = p3_lab,
                     do_statistics = F,
                     factor_colors = factor_colors
)
p0 = p0 %>% wrap_plots(nrow = 1)
p0
ggsave('fig/EFig10/kinematics.pdf', p0, width = 16, height = 8, units = 'cm', useDingbats = F)

sc = readRDS('data/regen/seurat/regen_final.rds')
sc@meta.data %<>%
    mutate(
        replicate = paste0(label, '_', slide)
    )

expr = GetAssayData(sc, slot='counts')
n_genes = colSums(expr > 0)
n_umis = colSums(expr)
pct_mito = PercentageFeatureSet(sc, pattern = "^mt-") %>% pull(nCount_RNA)

# extract metadata
meta = sc@meta.data %>%
    mutate(n_genes = n_genes, 
           n_umis = n_umis,
           pct_mito = pct_mito) %>%
    mutate(
        label = fct_recode(
            label, 
            'Young' = 'young', 
            'Old' = 'old',
            'Treated' = 'treated'
            ) %>% 
            fct_relevel('Young', 'Old', 'Treated'),
        replicate = fct_recode(
            replicate, 
            'old_1' = 'old_4', 
            'treated_2' = 'treated_4',
            'young_2' = 'young_4'
        ) %>% 
            fct_relevel(
                'young_1',
                'young_2',
                'young_3',
                'old_1',
                'old_2',
                'old_3',
                'treated_1',
                'treated_2'
            )
    )

mean_umis_per_barcode = meta %>%
    group_by(replicate, label) %>% 
    summarise(mean = mean(n_umis), median = median(n_umis),
              sd = sd(n_umis), sem = sd / sqrt(n())) %>% 
    ungroup() %>% 
    arrange(mean)
rep_order = mean_umis_per_barcode$replicate
p1a = mean_umis_per_barcode %>%
    arrange(mean) %>% 
    mutate(replicate = factor(replicate, levels=rep_order)) %>%
    ggplot(aes(x = mean, moe = sem, y = replicate, fill = label)) +
    stat_confidence_density(height = 0.8, confidence = 0.68, show.legend = TRUE) +
    geom_errorbarh(aes(xmin = mean, xmax = mean), size = 0.3, width = 0.8) +
    scale_x_continuous(expression('UMIs/\nbarcode'~(10^3)), 
                       labels = function(x) x / 1e3, limits = c(0, 20.5e3),
                       expand = c(0, 0)) +
    scale_y_discrete('Sections') +
    boxed_theme() +
    theme(legend.position = 'none',
          legend.key.size = unit(0.4, 'lines'),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_blank())
p1a

# number of genes
mean_genes_per_barcode = meta %>%
    group_by(replicate, label) %>% 
    summarise(mean = mean(n_genes), median = median(n_genes),
              sd = sd(n_genes), sem = sd / sqrt(n())) %>% 
    ungroup()
p1b = mean_genes_per_barcode %>%
    arrange(mean) %>% 
    mutate(replicate = factor(replicate, levels=rep_order)) %>%
    ggplot(aes(x = mean, moe = sem, y = replicate, fill = label)) +
    stat_confidence_density(height = 0.8, confidence = 0.68, show.legend = TRUE) +
    geom_errorbarh(aes(xmin = mean, xmax = mean), size = 0.3, width = 0.8) +
    scale_x_continuous(expression('Genes/\nbarcode'~(10^3)), 
                       labels = function(x) x / 1e3, limits = c(0, 5.5e3),
                       expand = c(0, 0)) +
    scale_y_discrete('Sections') +
    boxed_theme() +
    theme(legend.position = 'top',
          legend.key.size = unit(0.4, 'lines'),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
p1b

# percent mito
pct_mito_per_barcode = meta %>%
    group_by(replicate, label) %>% 
    summarise(mean = mean(pct_mito), median = median(pct_mito),
              sd = sd(pct_mito), sem = sd / sqrt(n())) %>% 
    ungroup()
p1c = pct_mito_per_barcode %>%
    arrange(mean) %>% 
    mutate(replicate = factor(replicate, levels=rep_order)) %>% 
    ggplot(aes(x = mean, moe = sem, y = replicate, fill = label)) +
    stat_confidence_density(height = 0.8, confidence = 0.68, show.legend = TRUE) +
    geom_errorbarh(aes(xmin = mean, xmax = mean), size = 0.3, width = 0.8) +
    scale_x_continuous(expression('% mitochondrial counts/\nbarcode'), 
                       limits = c(0, 21),
                       expand = c(0, 0)) +
    scale_y_discrete('Sections') +
    boxed_theme() +
    theme(legend.position = 'none',
          legend.key.size = unit(0.4, 'lines'),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
p1c

# combine and save together
p1 = p1a | p1b | p1c
p1
ggsave("fig/EFig10/per-section-metrics.pdf", p1,
       width = 7, height = 3.5, units = "cm", useDingbats = FALSE)

# histogram: number of genes
median1 = median(meta$n_genes)
p2_1 = meta %>%
    ggplot(aes(x = n_genes)) +
    geom_histogram(bins = 30, alpha = 0.4, size = 0.3) + 
    scale_x_continuous('# of genes, thousands', labels = function(x) x / 1e3,
                       breaks = seq(0, 10, 2) * 1e3) +
    scale_y_continuous('Barcodes', expand = c(0, 0), limits = c(0, 300),
                       breaks = seq(0, 500, 100)) +
    boxed_theme() +
    theme(legend.position = 'none') +
    geom_vline(aes(xintercept = median1), color = 'black', size = 0.3,
               linetype = 'dotted') +
    geom_label(data = data.frame(),
               aes(x = median1, y = Inf,
                   label = paste(format(median1, big.mark = ','), 'genes')),
               hjust = 0, vjust = 1, size = 1.5, fill = NA, label.size = NA,
               label.padding = unit(0.45, 'lines'))
p2_1
# ggsave("fig/EFig10/n-genes-histogram.pdf", p2_1, width = 4, height = 3, units = "cm", useDingbats = FALSE)
# 
# histogram: number of UMIs

median2 = median(meta$n_umis)
p2_2 = meta %>%
    ggplot(aes(x = n_umis)) +
    geom_histogram(bins = 30, alpha = 0.4, size = 0.3) + 
    scale_x_continuous('# of UMIs, thousands', labels = function(x) x / 1e3, limits=c(0, 30*1e3),
                       breaks = seq(0, 30, 10) * 1e3) +
    scale_y_continuous('Barcodes', expand = c(0, 0), limits = c(0, 300),
                       breaks = seq(0, 500, 100)) +
    boxed_theme() +
    theme(legend.position = 'none') +
    geom_vline(aes(xintercept = median2), color = 'black', size = 0.3,
               linetype = 'dotted') +
    geom_label(data = data.frame(),
               aes(x = median2, y = Inf,
                   label = paste(format(median2, big.mark = ','), 'UMIs')),
               hjust = 0, vjust = 1, size = 1.5, fill = NA, label.size = NA,
               label.padding = unit(0.45, 'lines'))
p2_2
# ggsave("fig/EFig10/n-UMIs-histogram.pdf", p2_2, width = 4, height = 3, units = "cm", useDingbats = FALSE)

# percent mitochondrial
median3 = round(median(meta$pct_mito),2)
p2_3 = meta %>%
    ggplot(aes(x = pct_mito)) +
    geom_histogram(bins = 30, alpha = 0.4, size = 0.3) + 
    scale_x_continuous('% mitochondrial counts', limits=c(0,30),
                       breaks = seq(0, 30, 10)) +
    scale_y_continuous('Barcodes', expand = c(0, 0), limits = c(0, 300),
                       breaks = seq(0, 500, 100)) +
    boxed_theme() +
    theme(legend.position = 'none') +
    geom_vline(aes(xintercept = median3), color = 'black', size = 0.3,
               linetype = 'dotted') +
    geom_label(data = data.frame(),
               aes(x = median3, y = Inf,
                   label = paste(format(median3, big.mark = ','), '% mito')),
               hjust = 0, vjust = 1, size = 1.5, fill = NA, label.size = NA,
               label.padding = unit(0.45, 'lines'))
p2_3
# ggsave("fig/EFig10/pct-mito-histogram.pdf", p2_3, width = 4, height = 3, units = "cm", useDingbats = FALSE)

# correlation
meta %<>% 
    mutate(dens = get_density(n_umis, n_genes))
range = range(meta$dens)
lims = c(range[1] + 0.1 * diff(range),
         range[2] - 0.1 * diff(range))
labs = c('min', 'max')

p2_4 = meta %>%
    arrange(dens) %>% 
    ggplot(aes(x = n_umis, y = n_genes, color = dens)) +
    rasterise(geom_point(size = 0.001, stroke = 0.25, shape = 20), dpi = 600) +
    scale_x_continuous('# of UMIs, thousands', labels = function(x) x / 1e3) +
    scale_y_continuous('# of genes, thousands', labels = function(x) x / 1e3,
                       breaks = seq(0, 10, 2) * 1e3) +
    # scale_color_paletteer_c("pals::linearl", name = 'Density', limits = range,
    #                        breaks = lims, labels = labs) +
    scale_color_gradientn(name = 'Density', colours = nr_heat_red_no_white,
                          limits = range, breaks = lims, labels = labs) +
    guides(color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
    boxed_theme() +
    theme(
        legend.position = 'bottom',
        legend.justification = 'right',
        legend.key.width = unit(0.15, 'lines'),
        legend.key.height = unit(0.15, 'lines')
        )
# p2_4
# ggsave("fig/EFig10/n-UMIs-vs-n-genes-density.pdf", p2_4, width = 4, height = 4.33, units = "cm", useDingbats = FALSE)

p2 = wrap_plots(p2_1,p2_2,p2_3,p2_4, ncol=4)
ggsave("fig/EFig10/summary-stats.pdf", p2, width = 11, height = 3.5, units = "cm", useDingbats = FALSE)

p3 = wrap_plots(p1, p2, widths=c(7,4))
# p3
ggsave("fig/EFig10/top-row.pdf", p3, width = 18, height = 12, units = "cm", useDingbats = FALSE)

expr = GetAssayData(sc, slot='counts') %>% NormalizeData()
meta = sc@meta.data %>%
    mutate(
        label = fct_recode(
            label, 
            'Young' = 'young', 
            'Old' = 'old',
            'Treated' = 'treated'
        ) %>% as.character()
    )

go_df = readRDS('data/metadata/go_names.rds') %>%
    mutate(go_sub = gsub('\\:', '-', go))

comparisons = list(
    'young-old'=list(
        go_suffix = 'young_old',
        factor_lvls = c('Young', 'Old'),
        genes = c(
            'Spp1',
            'Actb',
            'Tmsb4x',
            'Tyrobp',
            'C1qc',
            'Snap25',
            'Clu'
        ),
        gos = c(
            'neutrophil chemotaxis',
            'proteolysis',
            'astrocyte projection',
            'apoptotic cell clearance',
            'presynapse'
        )
    ),
    'treated-old'=list(
        go_suffix = 'treated_old',
        factor_lvls = c('Old', 'Treated'),
        genes = c(
            'Mobp',
            'Plekhb1',
            'Cd68',
            'Ctsb'
        ),
        gos = c(
            'apoptotic cell clearance',
            'immune response',
            'astrocyte development'
        )
    )
)

seurat_go_dir = 'data/regen/seurat_GO/'
for (comparison in names(comparisons)) {
    genes = comparisons[[comparison]]$genes
    meta0 = meta %>%
        filter(label %in% comparisons[[comparison]]$factor_lvls) %>%
        mutate(label = fct_relevel(label, comparisons[[comparison]]$factor_lvls))
    for (gene in genes) {
        if (gene == first(genes)) plot_list = list()
        meta0$expr = expr[gene, meta0$barcode]
        ## by condition
        conditions = levels(meta0$label)
        coords = map_dfr(conditions, ~ {
            condition = .x
            fit = loess(expr ~ x * y, data = filter(meta0, label == condition),
                        span = 0.02, degree = 1)
            coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5)) %>% 
                mutate(interp = predict(fit, newdata = .)) %>% 
                drop_na(interp) %>% 
                mutate(condition = condition)
        })
        # range = range(coords$interp)
        range = quantile(coords$interp, c(0.01, 0.99))
        brks = c(range[1] + 0.1 * diff(range),
                 range[2] - 0.1 * diff(range))
        labels = c('min', 'max')
        pal = nr_heat_red_no_white %>% tail(-5)
        p1a = coords %>% 
            mutate(
                interp = winsorize(interp, quantile(coords$interp, c(0.01, 0.99))),
                gene = gene
            ) %>%
            ggplot(aes(x = x, y = y, fill = interp)) +
            facet_grid(~condition) +
            geom_raster() +
            scale_y_reverse(expand = c(0, 0), breaks = c(-150, 150),
                            labels = c(expression('R'%->%""),
                                       expression(""%<-%'L')),
                            name = gene
            ) +
            scale_x_continuous(expand = c(0, 0), breaks = c(-400, 400),
                               labels = c(expression(""%<-%'Rostral'), 
                                          expression('Caudal'%->%""))
            ) +
            scale_fill_gradientn(name = 'Expression', colors = pal,
                                 limits = range, breaks = brks, labels = labels) +
            guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
            coord_fixed() +
            ggtitle(gene) +
            boxed_theme() +
            theme(
                plot.title = element_text(size = 6, margin = margin(0,0,-5,0)),
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
        plot_list[[length(plot_list)+1]] = p1a
    }
    
    go_sc = readRDS(paste0(seurat_go_dir, 'regen_final_', comparisons[[comparison]]$go_suffix, '.rds'))
    go_expr = GetAssayData(go_sc)
    gos = comparisons[[comparison]]$gos
    
    meta0 %<>% mutate(barcode = gsub('-', '_', barcode))
    
    for (go in gos) {
        go_sub = go_df %>%
            filter(go_name == !!go) %>%
            pull(go_sub)
        
        meta0$expr = go_expr[go_sub, meta0$barcode]
        ## by condition
        conditions = levels(meta0$label)
        coords = map_dfr(conditions, ~ {
            condition = .x
            fit = loess(expr ~ x * y, data = filter(meta0, label == condition),
                        span = 0.02, degree = 1)
            coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5)) %>% 
                mutate(interp = predict(fit, newdata = .)) %>% 
                drop_na(interp) %>% 
                mutate(condition = condition)
        })
        # range = range(coords$interp)
        range = quantile(coords$interp, c(0.01, 0.99))
        brks = c(range[1] + 0.1 * diff(range),
                 range[2] - 0.1 * diff(range))
        labels = c('min', 'max')
        pal = nr_heat_blue_spatial
        p1a = coords %>% 
            mutate(
                interp = winsorize(interp, quantile(coords$interp, c(0.01, 0.99))),
                gene = gene
            ) %>%
            ggplot(aes(x = x, y = y, fill = interp)) +
            facet_grid(~condition) +
            geom_raster() +
            scale_y_reverse(expand = c(0, 0), breaks = c(-150, 150),
                            labels = c(expression('R'%->%""),
                                       expression(""%<-%'L')),
                            name = gene
            ) +
            scale_x_continuous(expand = c(0, 0), breaks = c(-400, 400),
                               labels = c(expression(""%<-%'Rostral'), 
                                          expression('Caudal'%->%""))
            ) +
            scale_fill_gradientn(name = 'Expression', colors = pal,
                                 limits = range, breaks = brks, labels = labels) +
            guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
            coord_fixed() +
            ggtitle(go) +
            boxed_theme() +
            theme(
                plot.title = element_text(size = 6, margin = margin(0,0,-5,0)),
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
        plot_list[[length(plot_list)+1]] = p1a
    }
    
    p1 = wrap_plots(plot_list, ncol=2)
    p1
    
    ggsave(paste0('fig/EFig10/', comparison, '-genes-gos.pdf'), p1, width = 8, height = ifelse(comparison == 'treated-old', 7, 10), units='cm')
}

de_res1 = readRDS('data/regen/regen_final_treated_old-seed=42-nsub=10.rds')$de_feature_result
de_res2 = readRDS('data/regen/regen_final_young_old-seed=42-nsub=10.rds')$de_feature_result

combined_df = de_res1 %>%
    mutate(p_val = -log(p_val)) %>%
    dplyr::select(gene, p_val) %>%
    dplyr::rename(
        compar1_pval = p_val
    ) %>%
    left_join(
        de_res2 %>%
            mutate(p_val = -log(p_val)) %>%
            dplyr::select(gene, p_val) %>%
            dplyr::rename(
                compar2_pval = p_val
            )
    ) %>%
    filter(
        !is.na(compar1_pval),
        !is.na(compar2_pval)
    )
compar1_quan_range = quantile(combined_df$compar1_pval, c(0.01, 0.99))
compar2_quan_range = quantile(combined_df$compar2_pval, c(0.01, 0.99))

p5 = combined_df %>%
    ggplot(aes(x=compar1_pval, y=compar2_pval)) +
    rasterise(geom_point(size = 0.2, shape = 21, stroke = 0, alpha = 1, fill='black', color='black'), dpi = 600) +
    boxed_theme() +
    xlim(compar1_quan_range) +
    ylim(compar1_quan_range) +
    xlab(expression(-log[10](P)~"[treated vs. old]")) +
    ylab(expression(-log[10](P)~"[young vs. old]")) +
    geom_abline(slope=1, intercept=0, size=0.1, linetype = 'dashed') +
    stat_cor(size=1.5) +
    theme(
        aspect.ratio = 1
    )
p5
ggsave('fig/EFig10/p_val_cor.pdf', p5, width=4, height=3, units='cm')
