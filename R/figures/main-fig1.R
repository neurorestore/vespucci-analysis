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

###############################################################################-
## Simulation setup ####
###############################################################################-

# read data
sc = readRDS("data/simulations/objects/input=circle_2x-ct=none-de_prob=0.2-de_size=2-ct_prob=0-ct_size=0-overlap=0-reps=3-rep_de_jitter=1-rep_depth_jitter=0.3-seed=0.rds")
meta = sc@meta.data

# extract data
dat0 = meta %>% 
    dplyr::select(barcode, x, y, label, label_sim, replicate)

# plot label
label_pal = c('#F4F27E', '#372367')
p1_1 = dat0 %>% 
    ggplot(aes(x = y, y = x, fill = label)) +
    # facet_grid(~ 'Condition') +
    rasterise(geom_point(shape = 21, stroke = NA, size = 0.3), dpi = 300) +
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

truth_pal = pals::kovesi.linear_blue_95_50_c20(100)
brks = c(1, 2)
labels = c('min', 'max')
p1_2 = dat0 %>%
    mutate(
        label_sim = ifelse(label_sim == 'label1', 'background', label_sim),
        label_sim = factor(label_sim, levels=c('background', 'label2')),
        label_sim = as.integer(label_sim)
    ) %>%
    # arrange(label_sim) %>%
    ggplot(aes(x = y, y = x, fill = label_sim, color = label_sim)) +
    ggtitle('Ground truth') +
    ggrastr::rasterise(
        geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_gradientn(colours = truth_pal, name = 'Perturbation   ',
                          breaks = brks, labels = labels) +
    scale_fill_gradientn(colours = truth_pal, name = 'Perturbation   ', 
                         breaks = brks, labels = labels) +
    guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE,
                                 title.position = "left"),
           color = guide_colorbar(frame.colour = 'black', ticks = FALSE, 
                                  title.position = "left")) +
    coord_fixed() +
    boxed_theme(size_sm = 5, size_lg = 5) +
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
p1_2

# p1 = p1_1 | p1_2
p1 = p1_2
ggsave(paste0("fig/final/Fig1/simulation-setup.pdf"), p1, width = 3, height = 4, units = "cm", useDingbats = FALSE)

#############################################################################-
## Vespucci AUC ####
#############################################################################-

auc_res = readRDS('data/simulations/vespucci/input=circle_2x-ct=none-de_prob=0.2-de_size=2-ct_prob=0-ct_size=0-overlap=0-reps=3-rep_de_jitter=1-rep_depth_jitter=0.3-seed=0-ves_seed=42-max_cells=100-lo=eu_ma.rds')$aucs

dat0 %<>%
    left_join(auc_res)

# interpolate in 2D
fit = loess(auc ~ x * y, data = dat0, span = 0.015)
dat0$auc_fit = predict(fit, dat0)

range = range(dat0$auc_fit)
brks = c(range[1] + 0.1 * diff(range),
         range[2] - 0.1 * diff(range))
labels = format(range, digits = 2)
labels = c(paste0(labels[1], ' '),
           paste0(' ', labels[2]))

auc_pal = cet_pal(100, name = 'l19') %>% rev()

p2 = dat0 %>%
    # arrange(auc_fit) %>%
    ggplot(aes(x = y, y = x, fill = auc_fit, color = auc_fit)) +
    ggtitle('Vespucci') +
    ggrastr::rasterise(
        geom_point(size = 0.3, shape = 21, stroke = 0), dpi = 300
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
    coord_fixed() +
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
          legend.position = 'right',
          legend.justification = 'bottom',
          legend.key.width = unit(0.18, 'lines'),
          legend.key.height = unit(0.18, 'lines'),
          # plot.title = element_text(size = 5)
          plot.title = element_blank()
          )
# p2
ggsave(paste0("fig/final/Fig1/simulation-AUC.pdf"), p2, width = 4, height = 4, units = "cm", useDingbats = FALSE)

# read region accuracies ####
reg = readRDS('data/simulations/summaries/spatial_acc.rds') %>%
    type_convert() %>%
    filter(
        input == 'circle_2x',
        type == 'AUPRC'
    ) %>%
    mutate(
        method = ifelse(method == 'vespucci', 'Meta-learning', 'Exhaustive search')
    )

# read timing results ####
# summarise Magellan
magellan_files = list.files('data/simulations/summaries/magellan', 
                            full.names = TRUE) %>% 
    extract(grepl('circle_2x', .))
magellan_time = map_df(magellan_files, ~ {
    readRDS(.x) %>%
        group_by(input, seed) %>%
        summarise(
            time = sum(time / 3600)
        ) %>%
        mutate(units = 'hours', method = 'Exhaustive search')
})
# magellan runs on 4 threads
magellan_time$time %<>% `*`(4)

# summarise Vespucci
vespucci_files = list.files('data/simulations/vespucci',
                            full.names = TRUE) %>% 
    extract(grepl('circle_2x', .))
vespucci_files_df = map_df(vespucci_files, ~ {
    ori_filename = .x
    filename = gsub('\\.rds|\\.csv', '', basename(ori_filename))
    names = c('ori_filename')
    values = c(ori_filename)
    
    for (item in strsplit(filename, '-')[[1]]){
        temp = strsplit(item, '=')[[1]]
        if (length(temp) > 1){
            names = c(names, temp[1])
            values = c(values, temp[2])
        }
    }
    names(values) = names
    values
})
vespucci_files = vespucci_files_df %>%
    filter(
        (
            (de_prob == 0.5 & input %in% c('stripes', 'circle_overlap', 'flag')) |
                (de_prob == 0.2 & input %in% c('gradient', 'radial', 'circle_2x'))
        ),
        max_cells == 100
    ) %>% pull(ori_filename)
vespucci_time = map_df(vespucci_files, function(x){
    tmp = readRDS(x)
    input = gsub('input=', '', basename(x))
    input = gsub('-.*', '', input)
    seed = gsub('.*-seed=', '', basename(x))
    seed = gsub('-.*', '', seed)
    full_time = 
        sum(tmp$time_tracking$global$time) + 
        sum(tmp$time_tracking$auc$time)/60 * 8 + # vespucci runs on 8 threads
        sum(tmp$time_tracking$cor_convergence$model_time) +
        sum(tmp$time_tracking$cor_convergence$predict_time) 
    full_time = full_time/60
    return(data.frame(
        input = input,
        seed = seed,
        time = as.numeric(full_time),
        units = 'hours',
        method = 'Meta-learning'
    ))
})
time = rbind(
    magellan_time %>% type_convert(),
    vespucci_time %>% type_convert()
)

# read AUPR results #### 
# ves_res = readRDS('data/simulations/summaries/de_results/all_vespucci_de_auroc_summary.rds') %>%
#     type_convert() %>%
#     dplyr::select(-input_file, -lo) %>%
#     dplyr::rename(iter = ves_seed) %>%
#     filter(
#         input == 'circle_2x',
#         iter == 42,
#         de_model == 'nebula_nbgmm',
#         # de_prob == 0.5,
#         max_cells == 100
#     ) %>%
#     dplyr::select(input, seed, de_method, auprc_integral) %>%
#     mutate(input='circle')
# ves_res = readRDS('data/simulations/summaries/de_results/all_vespucci_de_auroc_summary_new.rds') %>% 
ves_res = readRDS('data/simulations/summaries/de_results/all_vespucci_de_auroc_summary.rds') %>%
    filter(input == 'circle', max_cells == 100, ves_seed == 42) %>% 
    dplyr::select(input, seed, auprc_integral) %>% mutate(de_method = 'vespucci')

splatter_de_stats = readRDS('data/simulations/summaries/de_results/splatter_de_results.rds') %>%
    dplyr::select(-ori_filename) %>%
    mutate(
        de_method = 'DE only',
        iter = 0,
        prop = -1,
        nsub = -1
    ) %>%
    filter(
        input == 'circle_2x'
    )
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
        input == 'circle'
    ) %>%
    dplyr::select(-ngenes, -ori_gene_size, -p_value_treatment)
other_de_stats %>% dplyr::select(input, de_method) %>% table()
# splatter_de_stats = splatter_de_stats[colnames(ves_res)]
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
    # splatter_de_stats,
    ves_res,
    other_de_stats
) %>% dplyr::rename(val = auprc_integral) %>% filter(de_method != 'somde')
DE %>% dplyr::select(seed, de_method) %>% table()
DE %<>% left_join(color_set, by = 'de_method')

stopifnot(nrow(time) == 20)
labs = time %>%
    group_by(input, method) %>%
    summarise(
        time = max(time),
        # label = paste0(round(median(time), 1), ' h')
        label = paste0(round(mean(time), 1), ' h')
    )

# paired t-test
delta = time %>% 
    group_by(seed) %>% 
    summarise(delta = time[method == 'Meta-learning'] - 
                  time[method == 'Exhaustive search']) %>% 
    ungroup()
stopifnot(nrow(delta) == 10)
pval = t.test(delta$delta)$p.value %>% 
    format(format = 'f', digits = 2) %>% 
    paste0('p = ', .)
pval_df = data.frame(label = pval)

range = boxplot(time ~ method, dat = time)$stats %>% range
og_pal = c('Exhaustive search' = nr_base_5[2],
           'Meta-learning' = nr_base_5[1])
p3_1 = time %>%
    ggplot(aes(x = method, y = time)) +
    # facet_grid(~ input) +
    geom_label(data = pval_df, aes(x = 1.5, y = -Inf, label = '____________________________'),
               size = 1.75, vjust = 0, label.size = NA, color = NA, fill = 'grey96',
               label.padding = unit(0.2, 'lines')) +
    geom_label(data = pval_df, aes(x = 1.5, y = -Inf, label = label),
               size = 1.5, vjust = 0, label.size = NA, fill = NA,
               label.padding = unit(0.45, 'lines')) +
    geom_boxplot(aes(color = method, fill = method), 
                 alpha = 0.4, width = 0.6, size = 0.35, outlier.shape = NA) +
    # geom_jitter(shape = 21, size = 0.4, stroke = 0.25, height = 0, 
    #             width = 0.15) +
    geom_text(dat = labs, aes(y = time, x = method, label = label),
              color = 'black',
              size = 1.5,
              vjust=-1) +
    boxed_theme() +
    scale_color_manual(values = og_pal) +
    scale_fill_manual(values = og_pal) +
    scale_x_discrete(labels = ~ gsub(" ", '\n', .)) +
    scale_y_continuous('Runtime, hours', breaks = pretty_breaks(4), 
                       expand = expansion(c(0.25, 0.25))) +
    coord_cartesian(ylim = range) +
    theme(
        aspect.ratio = 1.7,
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, lineheight = 0.8),
        legend.position = 'none',
        legend.justification = 'bottom',
        legend.title = element_blank()) 
p3_1

#############################################################################-
## Vespucci vs. Magellan, region AUC ####
#############################################################################-

stopifnot(nrow(reg) == 20)
labs = reg %>%
    group_by(method, input) %>%
    summarise(
        val = max(val),
        # label = round(median(val), 3)
        label = round(mean(val), 3)
    )

# paired t-test
delta = reg %>% 
    group_by(seed) %>% 
    summarise(delta = val[method == 'Meta-learning'] - 
                  val[method == 'Exhaustive search']) %>% 
    ungroup()
stopifnot(nrow(delta) == 10)
pval = t.test(delta$delta)$p.value %>% 
    format(format = 'f', digits = 2, scientific=T) %>% 
    paste0('p = ', .)
pval_df = data.frame(label = pval)

range = boxplot(val ~ method, dat = reg)$stats %>% range
og_pal = c('Exhaustive search' = nr_base_5[2],
           'Meta-learning' = nr_base_5[1])
p3_2 = reg %>%
    ggplot(aes(x = method, y = val)) +
    # facet_grid(~ input) +
    geom_label(data = pval_df, aes(x = 1.5, y = -Inf, label = '____________________________'),
               size = 1.75, vjust = 0, label.size = NA, color = NA, fill = 'grey96',
               label.padding = unit(0.2, 'lines')) +
    geom_label(data = pval_df, aes(x = 1.5, y = -Inf, label = label),
               size = 1.5, vjust = 0, label.size = NA, fill = NA,
               label.padding = unit(0.45, 'lines')) +
    geom_boxplot(aes(color = method, fill = method), 
                 alpha = 0.4, width = 0.6, size = 0.35, outlier.shape = NA) +
    # geom_jitter(shape = 21, size = 0.4, stroke = 0.25, height = 0, 
    #             width = 0.15) +
    geom_text(dat = labs, aes(y = val, x = method, label = label),
              color = 'black',
              size = 1.5,
              vjust=-1) +
    boxed_theme() +
    scale_color_manual(values = og_pal) +
    scale_fill_manual(values = og_pal) +
    scale_x_discrete(labels = ~ gsub(" ", '\n', .)) +
    scale_y_continuous('AUPRC', breaks = pretty_breaks(4), 
                       expand = expansion(c(0.25, 0.25))) +
    coord_cartesian(ylim = range) +
    theme(
        aspect.ratio = 1.7,
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, lineheight = 0.8),
        legend.position = 'none',
        legend.justification = 'bottom',
        legend.title = element_blank()) 
p3_2
p3 = p3_1 | p3_2
ggsave(paste0("fig/final/Fig1/circle-row.pdf"), p3, width = 6, height = 5, units = "cm", useDingbats = FALSE)

# set color palette
# pal = nr_base_11_light[1:11]
pal = pals::kelly(17)
names(pal) = unique(DE$color[!DE$color %in% c('Magellan', 'Vespucci')])
pal['Vespucci'] = nr_base_5[1]
alpha_pal = c('max'=1, 'min'=0.3)

# add OOT methods
oot_df = data.frame(x_name = c('SPARK', 'SPADE', 'trendsceek', 'GPCounts', 'BOOST-GP', 'BOOSTMI', 'scGCO', 'HSIC'), input = 'circle') %>% 
    mutate(color = x_name) %>% 
    filter(!x_name %in% DE$x_name)
DE %<>% bind_rows(oot_df)

# ensure all 10 simulations
for (de_method in unique(DE$de_method) %>% na.omit()) {
    print(de_method)
    rows = DE %>%
        filter(de_method == !!de_method) %>%
        nrow()
    stopifnot(rows == 10)
}

labs = DE %>%
    group_by(x_name, input, color) %>%
    summarize(
        # stats_val = median(val),
        stats_val = mean(val),
        val = max(val)
    ) %>%
    ungroup() %>%
    mutate(
        text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 2)),
        text_y = ifelse(val < 0, 0, val)
    ) %>% 
    replace_na(list(text_val = 'OOT', text_y = -Inf))

med = function(x) stats::median(x) %>% replace(is.na(.), 0)
p4 = DE %>%
    ggplot(aes(x = reorder(x_name, val, med), y = val,
               fill = color, color = color)) +
    geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
    geom_label(dat = labs,
               aes(label = text_val, y = text_y), 
               color = ifelse(labs$text_val == 'OOT', 'grey', 'black'),
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0, vjust = 0.5,
               show.legend = FALSE) +
    coord_flip() +
    scale_y_continuous('AUPR', breaks = pretty_breaks(), limits = c(0.18, 1),
                       expand = expansion(c(0.03, 0.125))) +
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
p4
ggsave(paste0("fig/final/Fig1/circle-AUPR-boxplot.pdf"), p4, width = 8, height = 7.3, units = "cm", useDingbats = FALSE)

#############################################################################-
## deltha heatmap ####
#############################################################################-

# test differences ####
pairs = tidyr::crossing(method1 = unique(DE$de_method),
                        method2 = unique(DE$de_method),
                        input = unique(DE$input)) %>% 
    filter(method1 != method2)
library(lawstat)
library(nparcomp)
tests = pmap_dfr(pairs, function(...) {
    current = tibble(...)
    print(current)
    
    vec1 = filter(DE, de_method == current$method1, input == current$input) %>% 
        arrange(input, seed) %>% pull(val)
    vec2 = filter(DE, de_method == current$method2, input == current$input) %>% 
        arrange(input, seed) %>% pull(val)
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
    vec1 = filter(DE, de_method == current$method1, input == current$input) %>% 
        arrange(seed) %>% pull(val)
    vec2 = filter(DE, de_method == current$method2, input == current$input) %>% 
        arrange(seed) %>% pull(val)
    median = median(vec2 - vec1)
    mean = mean(vec2 - vec1)
    mutate(current, delta_median = median, delta_mean = mean)
})

lvls = DE %>% 
    drop_na(val) %$%
    reorder(de_method, val, mean) %>% 
    levels()
range = range(deltas$delta_mean)
brks = c(range[1] + 0.1 * diff(range),
         range[2] - 0.1 * diff(range))
xlab = with(color_set, setNames(x_name, de_method))
labels = tests %>% 
    # filter(test == 'paired brunner-munzel') %>% 
    filter(test == 'paired t-test') %>% 
    mutate(method1 = factor(method1, levels = lvls),
           method2 = factor(method2, levels = lvls),
           lab = ifelse(pval < 0.001, '***',
                        ifelse(pval < 0.01, '**',
                               ifelse(pval < 0.05, '*', ''))))
p5 = deltas %>% 
    mutate(method1 = factor(method1, levels = lvls),
           method2 = factor(method2, levels = lvls)) %>% 
    ggplot(aes(x = method1, y = method2)) +
    # geom_tile(color = 'white', aes(fill = delta_median)) +
    geom_tile(color = 'white', aes(fill = delta_mean)) +
    geom_text(data = labels, size = 1.5, aes(label = lab), nudge_y = -0.15) +
    scale_x_discrete(expand = c(0, 0), labels = xlab) +
    scale_y_discrete(expand = c(0, 0), labels = xlab) +
    scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                           name = expression(Delta~AUPR),
                           breaks = brks,
                           labels = format(range, digits = 2)) +
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
p5
ggsave(paste0("fig/final/Fig1/circle-delta-heatmap.pdf"), p5, width = 9, height = 6.6, units = "cm", useDingbats = FALSE)