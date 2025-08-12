library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(sparseMatrixStats)

# setwd('C:/Users/teo/Documents/EPFL/projects/vespucci/')
setwd('~/git/vespucci/')
source('R/theme.R')

ves_de_summ = readRDS('/work/upcourtine/vespucci/rejected_review/find_usable_genes/vespucci_de_genes_summary.rds') %>% mutate(comparison = gsub('regen_final_', '', dataset)) %>%
    type_convert() %>%
    mutate(
        de_method_clean = case_when(
            de_method == 'sparkx' ~ 'SPARK-X',
            de_method == 'spacgn' ~ 'SpaGCN',
            de_method == 'cside' ~ 'C-SIDE',
            de_method == 'moransi' ~ 'Moransi test',
            de_method == 'wilcox' ~ 'Wilcoxon rank-sum test',
            de_method == 'nnsvg' ~ 'nnSVG',
            de_method == 'haystack' ~ 'SingleCellHaystack',
            de_method == 'spatialDE' ~ 'SpatialDE',
            de_method == 'spatialDE2' ~ 'SpatialDE2',
            de_method == 'heartsvg' ~ 'HEARTSVG',
            de_method == 'squidpy_permutation' ~ 'Squidpy (permutation test)',
            de_method == 'squidpy_normality' ~ 'Squidpy (normality ass.)',
            de_method == 'squidpy_normal_approx_permutation' ~ 'Squidpy (normal approx.)',
            de_method == 'scran' ~ 'scran',
            de_method == 'mast' ~ 'MAST',
            de_method == 'binSpect_kmeans' ~ 'binSpect (k-means)',
            de_method == 'binSpect_rank' ~ 'binSpect (rank)',
            de_method == 'dCor' ~ 'dCor',
            de_method == 'meringue' ~ 'MERINGUE',
            # de_method == 'somde' ~ 'SOMDE',
            de_method == 'spagft' ~ 'SpaGFT',
            de_method == 'spanve' ~ 'Spanve'
        )
    ) %>% 
    filter(!is.na(de_method_clean))

tmp_ves_de_summ = ves_de_summ %>% filter(comparison == 'young_old')
de_res = readRDS('/work/upcourtine/vespucci/real_data/DE_summaries/vespucci/regen_final_young_old-seed=42-nsub=10-de=nebula_nbgmm.rds') %>%
    mutate(pval_rank = rank(p_val_adj, ties='first')) %>% filter(gene %in% tmp_ves_de_summ$gene)
sc = readRDS('/work/upcourtine/vespucci/real_data/seurat/regen_final_young_old.rds')
meta = sc@meta.data %>% mutate(label = factor(str_to_title(label), levels=c('Old', 'Young')))
expr = GetAssayData(sc, slot='counts') %>% NormalizeData()

selected_genes = c('Spp1', 'C1qc', 'Lpl', 'Ctsd', 'B2m', 'Stab1')

plot_list = list()
for (gene in selected_genes) {
    meta0 = meta
    tmp_de_res = de_res %>% filter(gene == !!gene)
    rank_no = tmp_de_res$pval_rank
    meta0$expr = expr[gene, meta0$barcode]
    conditions = levels(meta0$label)
    coords = map_dfr(conditions, ~ {
        condition = .x
        fit = loess(expr ~ x * y, data = filter(meta0, label == condition), span = 0.02, degree = 1)
        coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5)) %>%
            mutate(interp = predict(fit, newdata = .)) %>%
            drop_na(interp) %>%
            mutate(condition = condition)
    })
    
    # range = range(coords$interp)
    range = quantile(coords$interp, c(0.01, 0.99))
    brks = c(range[1] + 0.1 * diff(range), range[2] - 0.1 * diff(range))
    labels = c('min', 'max')
    pal = nr_heat_red_no_white %>% tail(-5)
    p1 = coords %>%
        mutate(
            interp = winsorize(interp, quantile(coords$interp, c(0.01, 0.99))),
            gene = gene
        ) %>%
        ggplot(aes(x = x, y = y, fill = interp)) +
        facet_wrap(~condition, nrow=1) +
        geom_raster() +
        scale_y_reverse(expand = c(0, 0), breaks = c(-150, 150), labels = c(expression('R'%->%""), expression(""%<-%'L')), name = gene) +
        scale_x_continuous(expand = c(0, 0), breaks = c(-400, 400), labels = c(expression(""%<-%'Rostral'), expression('Caudal'%->%""))
        ) +
        scale_fill_gradientn(name = 'Expression', colors = pal, limits = range, breaks = brks, labels = labels) +
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
        coord_fixed() +
        ggtitle(str_to_title(gene)) +
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
    p1
    plot_list[[length(plot_list)+1]] = p1
}
p_out = wrap_plots(plot_list,ncol=3)
ggsave('fig/final/EFig16/chosen_gene_2d_plots.pdf', p_out, width=12, height=4, units='cm')

plot_list = list()
for (gene in selected_genes) {
    tmp_de_res = de_res %>% filter(gene == !!gene)
    rank_no = tmp_de_res$pval_rank
    tmp_ves_de_summ0 = tmp_ves_de_summ %>% filter(gene == !!gene) %>% dplyr::select(gene, de_method, de_method_clean, gene_rank) %>% rbind(data.frame(gene=gene, de_method='vespucci', de_method_clean = 'Vespucci', gene_rank=rank_no))
    tmp_ves_de_summ0 %<>% arrange(gene_rank)
    tmp_ves_de_summ0 %<>% mutate(de_method_clean = factor(de_method_clean, levels=rev(tmp_ves_de_summ0$de_method_clean)))
    pal = c(nr_heat_red[80]) %>% setNames(c('vespucci'))
    
    p2 = tmp_ves_de_summ0 %>%
        # filter(gene_rank < 1000) %>%
        mutate(gene_rank = ifelse(gene_rank >= 1000, 1000, gene_rank)) %>%
        mutate(group = factor(ifelse(de_method == 'vespucci', 'Vespucci', 'Others'), levels = c('Vespucci', 'Others')), rank_no = rank(gene_rank)) %>% 
        # mutate(gene_rank = log10(gene_rank)) %>%
        ggplot(aes(x = rank_no, y = gene_rank, color=de_method, fill=de_method)) +
        geom_point(size = 0.6, shape = 21, stroke = 0.2, position=position_dodge(width = .5)) +
        scale_y_continuous('Gene rank', breaks=seq(0, 1000, 200), labels=c(seq(0, 800, 200), '> 1000')) +
        ggtitle(paste0(gene, '; ', rank_no)) +
        facet_wrap(~group, nrow=2, strip.position = "left", scales = 'free_y') +
        scale_fill_manual(values=pal, na.value='grey20') +
        scale_color_manual(values=pal, na.value='grey20') +
        boxed_theme() +
        coord_flip() +
        theme(
            aspect.ratio = 0.1,
            strip.text.y.left = element_text(angle = 0),
            # axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            # axis.text.x = element_text(angle=45, hjust=1),
            # axis.ticks.x = element_blank(),
            # axis.text.x = element_blank(),
            legend.position = 'none'
        )
    plot_list[[length(plot_list)+1]] = p2
}

p_out = wrap_plots(plot_list, nrow=2)
ggsave('fig/final/EFig16/chosen_gene_rank.pdf', p_out, width=18, height=7, units='cm')


# p2 = tmp_ves_de_summ0 %>%
#     mutate(gene_rank = ifelse(gene_rank > 1000, 1000, gene_rank)) %>%
#     ggplot(aes(x = de_method_clean, y = gene_rank, color=de_method, fill=de_method)) +
#     geom_hline(aes(yintercept = 0), size = 0.4, color = 'grey88') +
#     geom_segment(aes(xend = de_method_clean, yend = 0)) +
#     geom_point(shape = 21, stroke = 0.2, size = 0.9) +
#     # scale_y_continuous(expression(paste(log[2], 'rank'))) +
#     scale_y_continuous('Gene rank') +
#     scale_fill_manual(values=pal, na.value='grey80') +
#     scale_color_manual(values=pal, na.value='grey80') +
#     # boxed_theme() +
#     clean_theme() +
#     coord_flip() +
#     theme(
#         # strip.text = element_blank(),
#         # aspect.ratio = 0.3,
#         axis.text.x = element_text(angle=45, hjust=1),
#         # axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = 'none'
#     )