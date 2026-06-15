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
source('R/theme.R')

sc = readRDS('data/simulations/objects/input=false-de_prob=0.5-de_size=2-reps=3-rep_de_jitter=1-rep_depth_jitter=0.3-sp_genes=f-seed=0.rds')
meta = sc@meta.data 
gene_features = sc@assays$originalexp@meta.features

# label_pal = c('#363062', '#435585', '#818FB4', '#F5E8C7')[c(2, 4)]
label_pal = c('#FFF5C2', '#F4F27E', '#6DB9EF', '#3081D0')[c(4, 2)]
p1 = meta %>% 
    ggplot(aes(x = y, y = x, fill = label)) +
    # facet_grid(~ 'Condition') +
    rasterise(geom_point(shape = 21, stroke = NA, size = 0.3), dpi = 600) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = label_pal, name = 'Condition', 
                      labels = c('Control', 'Treatment')) +
    guides(fill = guide_legend(override.aes = list(size = 1.0))) +
    # guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
    coord_fixed() +
    boxed_theme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0, 'lines'),
          axis.ticks.length.y = unit(0, 'lines'),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(rep(3.5, 4)),
          legend.position = 'top',
          # legend.justification = 'bottom',
          legend.key.width = unit(0.25, 'lines'),
          legend.key.height = unit(0.35, 'lines'),
          plot.title = element_text(size = 5))
# p1

effect_pal = pals::kovesi.linear_blue_95_50_c20(100)
p2 = meta %>% 
    ggplot(aes(x = y, y = x, fill = 1)) +
    # facet_grid(~ 'Condition') +
    rasterise(geom_point(shape = 21, stroke = NA, size = 0.3), dpi = 600) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colors = effect_pal, 
                      limits=c(1, 100),
                      breaks =c(1, 100),
                      name = 'Perturbation   ', 
                      labels = c('min', 'max')) +
    # guides(fill = guide_legend(override.aes = list(size = 1.0))) +
    coord_fixed() +
    guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE,
                                 title.position = "left"),
           color = guide_colorbar(frame.colour = 'black', ticks = FALSE, 
                                  title.position = "left")) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    # ggtitle('Perturbation Effect') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0, 'lines'),
          axis.ticks.length.y = unit(0, 'lines'),
          legend.position = 'top',
          legend.justification = 'bottom',
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.18, 'lines')
    )
p2

rep_pal = brewer.pal(6, 'Pastel2')
p3 = meta %>% 
    ggplot(aes(x = y, y = x, fill = replicate)) +
    # facet_grid(~ 'Condition') +
    rasterise(geom_point(shape = 21, stroke = NA, size = 0.3), dpi = 600) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = rep_pal, name = 'Library', 
                      labels = c(1:3, 1:3)) +
    guides(fill = guide_legend(override.aes = list(size = 1.0), nrow = 1)) +
    # guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
    coord_fixed() +
    boxed_theme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0, 'lines'),
          axis.ticks.length.y = unit(0, 'lines'),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(rep(3.5, 4)),
          legend.position = 'top',
          # legend.justification = 'bottom',
          legend.key.width = unit(0.25, 'lines'),
          legend.key.height = unit(0.35, 'lines'),
          plot.title = element_text(size = 5))

p0 = wrap_plots(p1, p2, p3, nrow=1)
ggsave("fig/final/EFig3/false-simulation-ground-truth.pdf", p0, width = 10, height = 4, units = "cm", useDingbats = FALSE)

plot_df = map_df(unique(meta$replicate), function(x){
    message(x)
    x_colname = gsub('Rep', 'Perturbation', x)
    vals = gene_features[,x_colname]
    data.frame(
        replicate=x,
        gene = gene_features$Gene,
        val = log10(vals)
    )
})

gene_truth = data.frame(
    'gene' = gene_features$Gene,
    'truth' = rowSums(gene_features[,endsWith(colnames(gene_features), 'is_selected')]) > 1
) %>%
    arrange(truth)

plot_df$gene = factor(plot_df$gene, levels=gene_truth$gene)

full_range = range(plot_df$val)
labels = round(full_range, 1)
color_pal = pals::kovesi.diverging_bwr_55_98_c37(100)
p2 = plot_df %>%
    ggplot(aes(x=replicate, y=gene, fill=val)) +
    geom_tile() +
    boxed_theme(size_lg = 4, size_sm = 4) +
    coord_flip() +
    scale_fill_gradientn(colours = color_pal, name='Effect size', breaks=full_range, labels=labels) +
    scale_colour_gradientn(colours = color_pal, name='Effect size', breaks=full_range, labels=labels) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    ylab('Gene') +
    ggtitle('Ground truth effect size') +
    theme(
        aspect.ratio = 0.3,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.14, 'lines'),
        plot.title = element_text(size=5)
    )

ggsave('fig/final/EFig3/false-simulations-effect.pdf', p2, width=8, height=6, units='cm')


source('R/theme.R')
ves_res = readRDS('data/simulations/summaries/de_results/vespucci_de_auroc_false_summary.rds') %>%
    type_convert() %>%
    dplyr::select(-input_file, -lo) %>%
    dplyr::rename(iter = ves_seed) %>%
    filter(
        # input %in% c('circle_overlap' ,'stripes', 'flag'),
        iter == 42,
        de_model == 'nebula_nbgmm',
        # de_prob == 0.5,
        max_cells == 100
    ) %>%
    dplyr::select(input, seed, de_method, fp)
splatter_de_stats = readRDS('data/simulations/summaries/de_results/splatter_de_false_results.rds') %>%
    dplyr::select(-ori_filename) %>%
    mutate(
        de_method = 'DE only',
        iter = 0,
        prop = -1,
        nsub = -1
    )
other_de_stats = readRDS('data/simulations/summaries/de_results/other_methods_false_stats.rds') %>%
    dplyr::select(-cell_type) %>% # no cell type effect for now
    mutate(
        iter = 0,
        prop = -1,
        nsub = -1
    ) %>%
    filter(
        de_method != 'smash', # not sure what to do with this for now,
        p_value_treatment == 'filtered'
    ) %>%
    dplyr::select(-ngenes, -ori_gene_size, -p_value_treatment)

other_de_stats %>% dplyr::select(input, de_method) %>% table()

splatter_de_stats = splatter_de_stats[colnames(ves_res)]
other_de_stats = other_de_stats[colnames(ves_res)]
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
dat0 = rbind(
    # splatter_de_stats,
    ves_res,
    other_de_stats
) %>% dplyr::rename(val = fp) %>%
    mutate(input = 'false') %>% filter(!de_method %in% c('somde', 'hsic'))

dat0 %>% dplyr::select(seed, de_method) %>% table()
dat0 %<>%
    left_join(color_set, by='de_method')

# set color palette
pal = pals::kelly(17)
names(pal) = unique(dat0$color[!dat0$color %in% c('Magellan', 'Vespucci')])
pal['Vespucci'] = nr_base_5[1]
alpha_pal = c('max'=1, 'min'=0.3)

# add OOT methods
oot_df = data.frame(x_name = c('SPARK', 'SPADE', 'trendsceek', 'GPCounts', 'BOOST-GP', 'BOOSTMI', 'scGCO', 'HSIC'),
                    input = dat0$input %>% first) %>% 
    mutate(color = x_name) %>% 
    filter(!x_name %in% dat0$x_name)
dat0 %<>% bind_rows(oot_df)

# ensure all 10 simulations
for (de_method in unique(dat0$de_method) %>% na.omit()) {
    print(de_method)
    rows = dat0 %>%
        filter(de_method == !!de_method) %>%
        nrow()
    stopifnot(rows == 10)
}

labs = dat0 %>%
    group_by(x_name, input, color) %>%
    summarize(
        stats_val = median(val),
        val = max(val)
    ) %>%
    ungroup() %>%
    mutate(
        text_val = as.character(ifelse(stats_val < 0, 'NA', format(stats_val, digits = 2))),
        text_y = ifelse(val < 0, 0, val)
    ) %>% 
    tidyr::replace_na(list(text_val = 'OOT', text_y = 0))

# dat0$color %<>% factor(levels = names(pal))
# dat0$alpha = ifelse(dat0$color %in% c('Magellan', 'Vespucci'), 'max', 'min')

x_name_levels = labs %>% arrange(-stats_val) %>% pull(x_name)
dat0$x_name = factor(dat0$x_name, levels = x_name_levels)

p5 = dat0 %>%
    ggplot(aes(x = x_name, y = val,
               fill = color, color = color)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = text_val, y = text_y), color = ifelse(labs$text_val == 'OOT', 'grey50', 'black'),
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0, vjust = 0.5,
               show.legend = FALSE) +
    coord_flip() +
    scale_y_continuous('# of false discoveries', breaks = pretty_breaks(),
                       expand = expansion(c(0.03, 0.125)), limits = c(0, 16000)) +
    scale_color_manual('', values = pal) +
    scale_fill_manual('', values = pal) +
    boxed_theme() +
    theme(axis.title.y = element_blank(),
          aspect.ratio = 1.5,
          legend.position = 'none',
          legend.justification = 'bottom',
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.15, 'lines')
    )
p5
ggsave(paste0("fig/final/EFig3/false-discoveries-boxplot.pdf"), p5,
       width = 8, height = 6.3, units = "cm", useDingbats = FALSE)

#############################################################################-
## f. delta-AUPR heatmap ####
#############################################################################-
pairs = tidyr::crossing(method1 = unique(dat0$de_method),
                        method2 = unique(dat0$de_method)) %>% 
    filter(method1 != method2)
library(lawstat)
library(nparcomp)
tests = pmap_dfr(pairs, function(...) {
    current = tibble(...)
    print(current)
    
    vec1 = filter(dat0, de_method == current$method1) %>% 
        arrange(seed) %>% pull(val)
    vec2 = filter(dat0, de_method == current$method2) %>% 
        arrange(seed) %>% pull(val)
    t = t.test(vec1, vec2)$p.value
    pt = t.test(vec1, vec2, paired = TRUE)$p.value
    w = wilcox.test(vec1, vec2)$p.value
    pw = wilcox.test(vec1, vec2, paired = TRUE)$p.value
    bm = brunner.munzel.test(vec1, vec2)$p.value
    pbm_df = data.frame(method = rep(c('method1', 'method2'), each = length(vec1)),
                        idx = rep(seq_along(vec1), 2),
                        value = c(vec1, vec2)) %>% 
        arrange(idx)
    pbm = npar.t.test.paired(value ~ method, pbm_df, nperm = 0)$Analysis['BM', 'p.value']
    data.frame(test = c('t-test', 'paired t-test', 'wilcox', 'paired wilcox',
                        'brunner-munzel', 'paired brunner-munzel'),
               pval = c(t, pt, w, pw, bm, pbm)) %>% 
        cbind(current, .)
})
# calculate delta-mean and median
deltas = pmap_dfr(pairs, function(...) {
    current = tibble(...)
    vec1 = filter(dat0, de_method == current$method1) %>% 
        arrange(seed) %>% pull(val)
    vec2 = filter(dat0, de_method == current$method2) %>% 
        arrange(seed) %>% pull(val)
    # median = median(vec2 - vec1)
    # mean = mean(vec2 - vec1)
    median = median(vec1 - vec2)
    mean = mean(vec1 - vec2)
    mutate(current, delta_median = median, delta_mean = mean)
})

delta0 = deltas %>%
    mutate(
        delta_median = -delta_median,
        delta_mean = -delta_mean
    )

lvls = dat0 %>% 
    drop_na(val) %>%
    group_by(de_method) %>%
    summarise(
        med = median(val)
        # med = mean(val)
    ) %>%
    ungroup() %>%
    arrange(-med, de_method) %>%
    pull(de_method)
# range = range(delta0$delta_median)
range = range(delta0$delta_mean)
brks = range
xlab = with(color_set, setNames(x_name, de_method))
labels = tests %>% 
    # filter(test == 'paired brunner-munzel') %>% 
    filter(test == 'paired t-test') %>% 
    mutate(method1 = factor(method1, levels = lvls),
           method2 = factor(method2, levels = lvls),
           lab = ifelse(pval < 0.001, '***',
                        ifelse(pval < 0.01, '**',
                               ifelse(pval < 0.05, '*', '')))
           ) %>%
    filter(as.integer(method2) < as.integer(method1))
delta0 %<>%
    mutate(
        method1 = factor(method1, levels = lvls),
        method2 = factor(method2, levels = lvls)
    ) %>%
    filter(
        ((as.integer(method2) < as.integer(method1)) | method2 == 'vespucci') 
        ) %>%
    mutate(
        delta_median = ifelse(method2 == 'vespucci', 0, delta_median),
        delta_mean = ifelse(method2 == 'vespucci', 0, delta_mean)
    )

p6 = delta0 %>% 
    ggplot(aes(x = method1, y = method2)) +
    geom_tile(color = 'white', aes(fill = delta_mean)) +
    # geom_tile(color = 'white', aes(fill = delta_median)) +
    geom_text(data = labels, size = 1.5, aes(label = lab), nudge_y = -0.15) +
    scale_x_discrete(expand = c(0, 0), labels = xlab) +
    scale_y_discrete(expand = c(0, 0), labels = xlab) +
    scale_fill_paletteer_c("pals::kovesi.linear_blue_95_50_c20",
                           name = expression(Delta~False~discoveries),
                           breaks = brks,
                           label = range
                           ) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    coord_fixed() +
    boxed_theme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          # axis.text.x = element_blank(),
          # axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.15, 'lines'),
          legend.position = 'right',
          legend.justification = 'bottom')
p6
ggsave(paste0("fig/final/EFig3/false-discoveries-delta.pdf"), p6,
       width = 9, height = 7, units = "cm", useDingbats = FALSE)
