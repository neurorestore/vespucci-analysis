setwd('C:/Users/teo/Documents/EPFL/projects/vespucci/')
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

###############################################################################-
## a. genes and GO expression plot ####
###############################################################################-
sc = readRDS('data/real_data/seurat/Koupourtidou2024.rds')
meta = sc@meta.data %>%
    mutate(barcode = gsub('-', '_', barcode))

expr = GetAssayData(sc, slot = 'counts')
colnames(expr) = gsub('-', '_', colnames(expr))

# normalize
expr %<>% NormalizeData()

# extract coordinates
dat0 = meta %>% dplyr::select(barcode, x, y, label)

# list genes to plot
genes_to_plot = c(
    'FTL1', 'C1QC', 'CTSZ', 'LY86', 'TIMP1', 'CCL12'
)
stopifnot(all(str_to_title(genes_to_plot) %in% rownames(expr)))
# genes_to_plot[!str_to_title(genes_to_plot) %in% rownames(expr)]

# iterate through genes

all_plots = list()
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
    
    pal = nr_heat_red_no_white %>% tail(-5)
    p = plot_df %>%
        mutate(
            label = ifelse(label == 'Intact', 'Intact', '3 dpi'),
            label = factor(label, levels=c('Intact', '3 dpi'))
        ) %>%
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
            strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0),
            aspect.ratio = 0.8,
            plot.title = element_text(size = 5, margin = margin(0,0,-10,0)),
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
            legend.justification = 'right'
        ) +
        facet_grid(~label) +
        coord_fixed()
    if (idx != length(genes_to_plot))
        p = p + theme(legend.position = 'none')
    p
    expr_plots[[idx]] = p
    p = p + theme(legend.position = 'none')
    all_plots[[length(all_plots)+1]] = p
}
p1 = wrap_plots(expr_plots, nrow = 2)
# ggsave('fig/final/EFig8/genes.pdf', p1, width = 10, height = 5,
#        units = 'cm', useDingbats = FALSE)


# load GO
sc = readRDS('data/real_data/seurat_GO/DE/Koupourtidou2024.rds')[[1]]
go_df = readRDS('data/metadata/go_names.rds')
meta = sc@meta.data %>%
    mutate(barcode = gsub('-', '_', barcode))
mat = GetAssayData(sc, slot='counts')
gos_to_plot = c(
    'positive regulation of microglial cell mediated cytotoxicity',
    'cysteine-type endopeptidase activity',
    'IgG binding'
)
gos = go_df %>% filter(go_name %in% gos_to_plot) %>%
    mutate(go_name = factor(go_name, levels=gos_to_plot)) %>%
    arrange(go_name) %>%
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
        mutate(
            label = ifelse(label == 'Intact', 'Intact', '3 dpi'),
            label = factor(label, levels=c('Intact', '3 dpi'))
        ) %>%
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
            strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0),
            aspect.ratio = 0.8,
            plot.title = element_text(size = 5, margin = margin(0,0,-10,0)),
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
            legend.justification = 'right'
        ) +
        facet_grid(~label) +
        coord_fixed()
    if (idx != length(gos_to_plot))
        p = p + theme(legend.position = 'none')
    go_plots[[idx]] = p
    p = p + theme(legend.position = 'none')
    all_plots[[length(all_plots)+1]] = p
}

p2 = wrap_plots(go_plots, nrow=1)
# ggsave('fig/final/EFig8/GO-modules.pdf', p2, width = 10, height = 3,
#        units = 'cm', useDingbats = FALSE)

p3 = wrap_plots(all_plots, nrow=3)
ggsave('fig/final/EFig8/genes-GO-modules.pdf', p3, width = 10, height = 6,
       units = 'cm', useDingbats = FALSE)


###############################################################################-
## b. lollipop plot: genes ####
###############################################################################-

# load results
dat0 = readRDS('data/real_data/DE_summaries/vespucci/Koupourtidou2024-seed=42-nsub=10-de=nebula_nbgmm.rds')
min_pval = min(dat0$p_val[dat0$p_val > 0])
dat0$p_val = vapply(dat0$p_val, function(x){max(min_pval, x)}, as.numeric(1))
dat0 %<>%
    mutate(gene = ifelse(gene == '1500015O10Rik', 'Ecrg4', gene)) %>%
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

p4 = dat0 %>%
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
p4
ggsave('fig/final/EFig8/lollipop-genes.pdf', p4, width = 8, height = 8, 
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## c. lollipop plot: GO modules ####
###############################################################################-

# load GO
go_df = readRDS("data/metadata/go_names.rds")

# load GLM results
dat0 = readRDS('data/real_data/DE_summaries/vespucci_GO/Koupourtidou2024-seed=42-nsub=10-de=glm.rds') %>% 
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
    left_join(go_df, by=c('gene'='go'))
dat0 %<>% 
    arrange(-logFC)
dat0$go_name = factor(dat0$go_name, levels=rev(as.character(dat0$go_name)))

# plot
p5 = dat0 %>%
    ggplot(aes(x = go_name, y = logFC)) +
    # facet_wrap(sign(logFC) ~ ., ncol = 1, scales = 'free') +
    geom_hline(aes(yintercept = 0), size = 0.4, color = 'grey88') +
    geom_segment(aes(xend = go_name, yend = 0), color = 'grey88') +
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
p5
ggsave('fig/final/EFig8/lollipop-GO.pdf', p4, width = 18, height = 8, 
       units = 'cm', useDingbats = FALSE)

ves_res = readRDS('data/real_data/DE_summaries/vespucci/Koupourtidou2024-seed=42-nsub=10-de=nebula_nbgmm.rds')
sparkx_res = readRDS('data/real_data/DE_summaries/others/Koupourtidou2024-de=sparkx.rds')

res_rank_check = ves_res %>%
    mutate(ves_rank = rank(p_val_adj)) %>%
    arrange(p_val_adj) %>%
    dplyr::select(gene, ves_rank, p_val_adj) %>%
    dplyr::rename(ves_pval_adj = p_val_adj) %>%
    left_join(
        sparkx_res %>%
            mutate(sparkx_rank = rank(p_val_adj)) %>%
            arrange(p_val_adj) %>%
            dplyr::select(gene, sparkx_rank, p_val_adj) %>%
            dplyr::rename(sparkx_pval_adj = p_val_adj)
    )

ves_top_genes = genes_to_plot = res_rank_check %>% filter(ves_rank <= 12) %>% arrange(ves_rank) %>% pull(gene)
sparkx_top_genes = genes_to_plot = res_rank_check %>% filter(sparkx_rank <= 12) %>% arrange(sparkx_rank)%>% pull(gene)

pal = nr_heat_red_no_white %>% tail(-5)
types = c('sparkx')
for (type in types) {
    if (type == 'vespucci') {
        genes_to_plot = ves_top_genes
    } else {
        genes_to_plot = sparkx_top_genes
    }
    plot_list = list()
    
    for (gene in genes_to_plot) {
        dat0$expr = expr[gene,meta$barcode]
        
        dat2 = data.frame()
        for (label in unique(meta$label)) {
            dat_tmp = dat0 %>%
                filter(label == !!label)
            fit = loess(expr ~ x * y, data = dat_tmp, span = 0.02, degree = 1)
            dat2 %<>% rbind(dat_tmp %>%
                                mutate(
                                    expr_fit = predict(fit, .)
                                ) %>%
                                filter(!is.na(expr_fit))
            )
        }
        
        range = range(dat2$expr_fit)
        brks = c(range[1] + 0.1 * diff(range),
                 range[2] - 0.1 * diff(range))
        labels = c('min', 'max')
        
        p2_1 = dat2 %>%
            mutate(
                label = ifelse(label == 'Intact', 'Intact', '3 dpi'),
                label = factor(label, levels=c('Intact', '3 dpi'))
            ) %>%
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
            ggtitle(str_to_upper(gene)) +
            theme(
                strip.background = element_blank(),
                strip.text.y.left = element_text(angle = 0),
                aspect.ratio = 0.8,
                plot.title = element_text(size = 5, margin = margin(0,0,-10,0)),
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
                legend.justification = 'right'
            ) +
            facet_grid(~label) +
            coord_fixed()
        if (gene != last(genes_to_plot)) {
            p2_1 = p2_1 + theme(legend.position = 'none')
        }
        plot_list[[length(plot_list)+1]] = p2_1
    }
    
    p2 = wrap_plots(plot_list, nrow=3)
    ggsave(paste0('fig/final/EFig8/', type, '-top-genes.pdf'), p2, width=14.5, height=9.5, units = 'cm')
}

# pal = pals::kovesi.diverging_bwr_55_98_c37(100)
ranking_df = res_rank_check %>% 
    # filter(ves_rank <= 100) %>%
    mutate(
        type = 'Vespucci',
        rank = ves_rank
    ) %>%
    dplyr::select(gene, rank, type)
ranking_df %<>% rbind(
    res_rank_check %>%
        filter(gene %in% ranking_df$gene) %>%
        mutate(
            type = 'SPARK-X',
            rank = sparkx_rank
        ) %>%
        dplyr::select(gene, rank, type)
)

p3 = ranking_df %>% 
    ggplot(aes(x=type, y=1, fill=rank)) +
    geom_bar(stat='identity') + 
    scale_fill_gradientn(
        colours = pal, 
        breaks = range(ranking_df$rank)) +
    boxed_theme() + 
    ggtitle('Gene rank') +
    scale_y_continuous(expand=c(0,0)) +
    ylab('Gene rank based by p-value') +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black', title='Rank')) +
    theme(
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(),
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.14, 'lines'),
        legend.text = element_text(size = 5),
        legend.position = 'right',
        legend.justification = 'bottom'
        # legend.title = element_blank()
    )
p3
ggsave('fig/final/EFig8/sparkx-vespucci-ranking.pdf', p3, width=4, height=5, units='cm')
