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
sc = readRDS("data/simulations/objects/input=circle-de_prob=0.2-de_size=2-reps=3-rep_de_jitter=1-rep_depth_jitter=0.3-seed=0.rds")
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

ves_res = readRDS('data/simulations/vespucci/input=circle-de_prob=0.2-de_size=2-reps=3-rep_de_jitter=1-rep_depth_jitter=0.3-seed=0-ves_seed=42.rds')

auc_res = ves_res$spatial_auc_result$aucs
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
reg = readRDS('data/simulations/spatial-acc-summary.rds') %>%
    type_convert() %>%
    filter(
        input == 'circle',
        type == 'AUPRC'
    ) %>%
    mutate(
        method = ifelse(method == 'vespucci', 'Meta-learning', 'Exhaustive search')
    )

# read timing results ####
time = readRDS('data/simulations/time-summary.rds') %>%
    mutate(
        method = ifelse(method == 'vespucci', 'Meta-learning', 'Exhaustive search')
    ) %>%
    filter(
        input == 'circle'
    )

# read AUPR results #### 
ves_de_res = readRDS('data/simulations/vespucci_de/input=circle-de_prob=0.2-de_size=2-reps=3-rep_de_jitter=1-rep_depth_jitter=0.3-seed=0-ves_seed=42.rds')
other_de_stats = readRDS('data/simulations/summaries/de-results.rds') %>%
    filter(
        input == 'circle'
    )

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
        'heartsvg',
        'squidpy_permutation',
        'squidpy_normality',
        'squidpy_normal_approx_permutation',
        'scran',
        'mast',
        'binSpect_kmeans',
        'binSpect_rank',
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
        'HEARTSVG',
        rep('SquidPy',3),
        rep('Giotto', 4),
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
        'HEARTSVG',
        'Squidpy (permutation test)',
        'Squidpy (normality ass.)',
        'Squidpy (normal approx.)',
        'scran',
        'MAST',
        'binSpect (k-means)',
        'binSpect (rank)',
        'Vespucci'
    )
)

DE = rbind(
    ves_de_res,
    other_de_stats
) %>% dplyr::rename(val = auprc_integral)
DE %>% dplyr::select(seed, de_method) %>% table()
DE %<>% left_join(color_set, by = 'de_method')

labs = time %>%
    group_by(input, method) %>%
    summarise(
        time = max(time),
        label = paste0(round(mean(time), 1), ' h')
    )

# paired t-test
delta = time %>% 
    group_by(seed) %>% 
    summarise(delta = time[method == 'Meta-learning'] - 
                  time[method == 'Exhaustive search']) %>% 
    ungroup()

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
pal = nr_base_11_light[1:10]
names(pal) = unique(DE$color[!DE$color %in% c('Vespucci')])
pal['Vespucci'] = nr_base_5[1]
alpha_pal = c('max'=1, 'min'=0.3)

# add OOT methods
oot_df = data.frame(x_name = c('SPARK', 'SPADE', 'trendsceek', 'GPCounts'),
                    input = 'circle_2x') %>% 
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
               aes(label = text_val, y = text_y), color = ifelse(text_val == 'OOT', 'grey', 'black'),
               label.padding = unit(0.35, 'lines'),
               label.size = NA, fill = NA,
               size = 1.75, hjust = 0, vjust = 0.5,
               show.legend = FALSE) +
    coord_flip() +
    scale_y_continuous('AUPR', breaks = pretty_breaks(), limits = c(0.18, 0.9),
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
ggsave(paste0("fig/final/Fig1/circle-AUPR-boxplot.pdf"), p4, width = 8, height = 6.3, units = "cm", useDingbats = FALSE)


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
ggsave(paste0("fig/final/EFig3/circle-delta-heatmap.pdf"), p5, width = 9, height = 6.6, units = "cm", useDingbats = FALSE)


#############################################################################-
# Calcagno registration
#############################################################################-
sc = readRDS('data/published_data/seurat/Calcagno2022.rds')
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
auc_res = readRDS('data/published_data/vespucci/Calcagno2022-seed=42-nsub=10.rds')$spatial_auc_result$aucs
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
sc = readRDS('data/published_data/seurat/Calcagno2022.rds')
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
sc = readRDS('data/published_data/seurat_GO/Calcagno2022.rds')
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

#############################################################################-
# Koupourtidou auc
#############################################################################-

sc = readRDS('data/published_data/seurat/Koupourtidou2024.rds')
meta = sc@meta.data %>%
    mutate(barcode = gsub('-', '_', barcode))

auc_res = readRDS('data/published_data/vespucci/Koupourtidou2024-seed=42-nsub=10.rds')$spatial_auc_result$aucs
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

p10 = meta %>%
    # arrange(auc_fit) %>%
    ggplot(aes(x = x, y = y, fill = auc_fit, color = auc_fit)) +
    ggtitle('Vespucci') +
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
    coord_fixed() +
    boxed_theme(size_lg = 5, size_sm = 5) +
    theme(
        # aspect.ratio = 0.8,
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
        # plot.title = element_text(size = 5)
        plot.title = element_blank()
    )
# p10
ggsave(paste0("fig/final/Fig1/koupourtidou_AUC.pdf"), p10, width = 2, height = 2.7, units = "cm", useDingbats = FALSE)

###############################################################################-
## Koupourtidou genes ####
###############################################################################-
expr = GetAssayData(sc, slot = 'counts')
colnames(expr) = gsub('-', '_', colnames(expr))

# normalize
expr %<>% NormalizeData()

# extract coordinates
dat0 = meta %>% dplyr::select(barcode, x, y, label)

# list genes to plot
genes_to_plot = c(
    'LCN2', 'PTGDS'
)

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
            legend.justification = 'bottom'
        ) +
        facet_grid(label~.,switch = "y") +
        coord_fixed()
    if (idx != 2)
        p = p + theme(legend.position = 'none')
    if (idx != 1)
        p = p + theme(strip.text.y.left = element_blank())
    p
    expr_plots[[idx]] = p
}
# p11 = wrap_plots(expr_plots, ncol = 1)
p11 = wrap_plots(expr_plots, nrow = 1)
ggsave('fig/final/Fig1/koupourtidou_genes.pdf', p11, width = 3.5, height = 3.5,
       units = 'cm', useDingbats = FALSE)

###############################################################################-
## Koupourtidou GO terms ####
###############################################################################-

# load GO
sc = readRDS('data/published_data/seurat_GO/Koupourtidou2024.rds')[[1]]
go_df = readRDS('data/metadata/go_names.rds')
meta = sc@meta.data %>%
    mutate(barcode = gsub('-', '_', barcode))
mat = GetAssayData(sc, slot='counts')
gos_to_plot = c(
    'complement activation, classical pathway',
    'cysteine-type endopeptidase activity'
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
            label = ifelse(label == 'Intact', 'Intact', '3 day\npost injury'),
            label = factor(label, levels=c('Intact', '3 day\npost injury'))
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
            legend.justification = 'bottom'
        ) +
        facet_grid(label~.,switch = "y") +
        coord_fixed()
    if (idx != 1)
        p = p + theme(legend.position = 'none')
    if (idx != 2)
        p = p + theme(strip.text.y.left = element_blank())
    p
    go_plots[[idx]] = p
}
p12 = go_plots[[1]]
ggsave('fig/final/Fig1/koupourtidou_GO-modules.pdf', p12, width = 2.4, height = 3,
       units = 'cm', useDingbats = FALSE)

combined_plots = list(expr_plots[[1]], expr_plots[[2]], go_plots[[1]])
p13 = wrap_plots(combined_plots, nrow=1)
ggsave('fig/final/Fig1/koupourtidou_genes-and-GO-modules.pdf', p13, width = 5, height = 3.4,
       units = 'cm', useDingbats = FALSE)
