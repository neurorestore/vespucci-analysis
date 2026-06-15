library(tidyverse)
library(magrittr)
library(ggplot2)
setwd('~/git/vespucci-analysis/')
source('R/theme.R')

# read AUPR results #### 
ves_res = readRDS('data/simulations/summaries/de_results/all_vespucci_de_auroc_summary.rds') %>%
    type_convert() %>%
    dplyr::select(-input_file, -lo) %>%
    dplyr::rename(iter = ves_seed) %>%
    mutate(input = ifelse(input == 'circle_2x', 'circle', input)) %>% 
    filter(
        iter == 42,
        !input %in% c('gradient', 'radial'),
        de_model == 'nebula_nbgmm',
        # de_prob == 0.5,
        max_cells == 100,
        (de_prob == 0.5 & input %in% c('stripes', 'circle_overlap', 'flag')) |
		(de_prob == 0.2 & input %in% c('circle'))
    ) %>%
    dplyr::select(input, seed, de_method, auprc_integral, acc, sensitivity, specificity)
other_de_stats = readRDS('data/simulations/summaries/de_results/other_methods_stats.rds') %>%
    dplyr::select(-cell_type) %>% # no cell type effect for now
    mutate(
        iter = 0,
        prop = -1,
        nsub = -1
    ) %>%
    filter(
        de_method != 'hsic',
        p_value_treatment == 'filtered',
    ) %>%
    dplyr::select(-ngenes, -ori_gene_size, -p_value_treatment)
other_de_stats %>% dplyr::select(input, de_method) %>% table()
other_de_stats = other_de_stats[colnames(ves_res)]
color_set = data.frame(
    de_method = c(
        'DE only',
        'sparkx',
        'spacgn',
        'cside',
        'moransi',
        'wilcox',
        'nnsvg',
        'spatialDE',
        'spatialDE2',
        'heartsvg',
        'squidpy_permutation',
        'squidpy_normality',
        'squidpy_normal_approx_permutation',
        'scran',
        'mast',
        'binSpect_kmeans',
        'binSpect_rank',
        'haystack',
        'dCor',
        'hsic',
        'meringue',
        'rv',
        'somde',
        'spagft',
        'spanve',
        'magellan',
        'vespucci'
    ),
    color = c(
        'DE only',
        'SPARK-X',
        'SpaGCN',
        'C-SIDE',
        rep('Seurat',2),
        'nnSVG',
        'SpatialDE',
        'SpatialDE2',
        'HEARTSVG',
        rep('SquidPy',3),
        rep('Giotto', 4),
        'SingleCellHayStack',
        'dCor',
        'Hsic',
        'MERINGUE',
        'RV',
        'SomDE',
        'Spagft',
        'Spanve',
        'Magellan',
        'Vespucci'
    ),
    x_name = c(
        'NBGMM',
        'SPARK-X',
        'SpaGCN',
        'C-SIDE',
        'Moransi test',
        'Wilcoxon rank-sum test',
        'nnSVG',
        'SpatialDE',
        'SpatialDE2',
        'HEARTSVG',
        'Squidpy (permutation test)',
        'Squidpy (normality ass.)',
        'Squidpy (normal approx.)',
        'scran',
        'MAST',
        'binSpect (k-means)',
        'binSpect (rank)',
        'SingleCellHayStack',
        'dCor',
        'Hsic',
        'MERINGUE',
        'RV',
        'SomDE',
        'Spagft',
        'Spanve',
        'Magellan',
        'Vespucci'
    )
)
DE = rbind(
    ves_res,
    other_de_stats
) %>% dplyr::rename(auprc = auprc_integral)
DE %>% dplyr::select(seed, de_method) %>% table()
DE %<>% left_join(color_set, by = 'de_method')

pal = pals::kelly(17)
names(pal) = unique(DE$color[!DE$color %in% c('Magellan', 'Vespucci')])
pal['Vespucci'] = nr_base_5[1]
alpha_pal = c('max'=1, 'min'=0.3)

# add OOT methods
oot_df = tidyr::crossing(x_name = c('SPARK', 'SPADE', 'trendsceek', 'GPCounts', 'BOOST-GP', 'BOOSTMI', 'scGCO', 'HSIC', 'SOMDE'),input = unique(DE$input)) %>% 
    mutate(color = x_name) %>% 
    filter(!x_name %in% DE$x_name)
DE %<>% bind_rows(oot_df)

DE %<>% filter(de_method != 'somde')

# ensure all 10 simulations
for (de_method in unique(DE$de_method) %>% na.omit()) {
    print(de_method)
    rows = DE %>%
        filter(de_method == !!de_method) %>%
        nrow()
    stopifnot(rows == 40)
}

grid_to_run = tidyr::crossing(
    metric = c('acc', 'sensitivity', 'specificity'),
    input = unique(DE$input)
) %>% mutate(
    input = factor(input, levels=c('circle', 'circle_overlap', 'stripes', 'flag'))
) %>% 
arrange(metric, input)

plot_list = list()
for (i in 1:nrow(grid_to_run)) {
    tmp_row = grid_to_run[i,]
    DE0 = DE %>% filter(input == tmp_row$input) %>% mutate(val = get(tmp_row$metric), 
    val = ifelse(is.na(val) & !x_name %in% oot_df$x_name, 0, val))
    metric_name = case_when(
        tmp_row$metric == 'auprc' ~ 'AUPR',
        tmp_row$metric == 'acc' ~ 'Accuracy',
        tmp_row$metric == 'sensitivity' ~ 'Sensitivity',
        tmp_row$metric == 'specificity' ~ 'Specificity'
    )
    sim_name = paste0('Simulation #', as.integer(tmp_row$input))
    labs = DE0 %>%
        group_by(x_name, input, color) %>%
        summarize(
            # stats_val = median(val),
            stats_val = mean(val),
            val = max(val)
        ) %>%
        ungroup() %>%
        mutate(
            text_val = ifelse(stats_val < 0, 'NA', format(round(stats_val, 2))),
            text_y = ifelse(val < 0, 0, val)
        ) %>% 
        replace_na(list(text_val = 'OOT', text_y = -Inf))

    med = function(x) stats::median(x) %>% replace(is.na(.), 0)
    p0 = DE0 %>%
        ggplot(aes(x = reorder(x_name, val, med), y = val, fill = color, color = color)) +
        geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
        geom_label(dat = labs,
                aes(label = text_val, y = text_y), 
                color = ifelse(labs$text_val == 'OOT', 'grey', 'black'),
                label.padding = unit(0.35, 'lines'),
                label.size = NA, fill = NA,
                size = 1.75, hjust = 0, vjust = 0.5,
                show.legend = FALSE) +
        coord_flip() +
        scale_y_continuous(metric_name, breaks = pretty_breaks(), limits = c(0, 1), expand = expansion(c(0.03, 0.125))) +
        scale_color_manual('', values = pal) +
        scale_fill_manual('', values = pal) +
        boxed_theme() +
        ggtitle(sim_name) + 
        theme(axis.title.y = element_blank(),
            aspect.ratio = 1.5,
            legend.position = 'none',
            legend.justification = 'bottom',
            legend.key.width = unit(0.18, 'lines'),
            legend.key.height = unit(0.15, 'lines')
        )
    if (tmp_row$metric == 'specificity') {
        p0 = p0 + scale_y_continuous(metric_name, breaks = pretty_breaks(), limits = c(0, 1), expand = expansion(c(0.03, 0.2)))
    }
    plot_list[[length(plot_list)+1]] = p0
}

out_p = wrap_plots(plot_list, ncol=3, byrow=F)
ggsave(paste0("fig/EFig3/additional-metrics-boxplot.pdf"), out_p, width = 19, height = 22, units = "cm", useDingbats = FALSE)
