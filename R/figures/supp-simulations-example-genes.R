setwd('~/git/vespucci-analysis/')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(cetcolor)
library(patchwork)
library(ggrastr)
source('R/theme.R')

# inputs = c('circle', 'circle_overlap', 'stripes', 'flag')
input = 'circle'
de_methods = c(
    'vespucci',
    sort(
        c('binSpect_kmeans', 'binSpect_rank', 'cside', 'dCor', 'haystack', 'heartsvg', 'mast', 'meringue', 'moransi', 'nnsvg', 'rv', 'scran', 'sparkx', 'wilcox', 'spacgn', 'spagft', 'spanve', 'spatialDE', 'spatialDE2', 'squidpy')
    )
)
de_method_names = data.frame(
    de_method_name = c(
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

normalize_data = TRUE
input_prefix = 'input=circle-seed=0'
# sc_obj_filename = paste0('data/simulations/objects/',input_prefix,'.rds')
sc_obj_filename = paste0('data/simulations/objects/',input_prefix,'.rds')
sc = readRDS(sc_obj_filename)
genes = sc@assays$originalexp@meta.features
meta = sc@meta.data
expr = GetAssayData(sc, slot = 'counts')

# optionally, normalize
if (normalize_data) {
    expr %<>% NormalizeData()
    norm_suffix = 'norm'
} else {
    norm_suffix = 'raw'
}

# genes_meta = readRDS('data/simulations/objects_meta/input=circle-seed=0.rds')$gene_meta 
genes_meta = readRDS('data/simulations/objects_meta/input=circle-seed=0.rds')$gene_meta 
truth = genes_meta%>% 
    dplyr::select(colnames(genes_meta)[grepl('_is_selected', colnames(genes_meta))]) %>% 
    rowSums() %>%
    data.frame() %>%
    set_colnames('truth') %>%
    mutate(truth = as.integer(truth > 0)) %>% 
    rownames_to_column('gene')

# start with Vespucci first
for (de_method in de_methods) {
    if (de_method == first(de_methods)) {
        set.seed(42)
        expr_plots = list()
    }
    print(de_method)
    if (de_method != 'vespucci') {
        if (!de_method %in% c('spacgn', 'spagft', 'spanve', 'spatialDE', 'spatialDE2', 'squidpy')) {
            data_dir = 'data/simulations/DE_summaries/others/'
            input_suffix = '.rds'
        } else {
            data_dir = 'data/simulations/DE_summaries/others_python/'
            input_suffix = '.csv'
        }
        de_filename = paste0(data_dir, input_prefix, '-de=', de_method, input_suffix)
        
        if (endsWith(de_filename, '.csv')) {
            tmp = read.csv(de_filename)
        } else if (endsWith(de_filename, '.rds')) {
            tmp = readRDS(de_filename)
        }
    } else {
        tmp = readRDS('data/simulations/vespucci/input=circle-seed=0-ves_seed=42-max_cells=100.rds')$de_feature_result %>% dplyr::rename(gene = feature)
    }
    
    # manual fixes
    if (de_method == 'squidpy') {
        tmp %<>% mutate(gene = Gene)
    } else if (de_method == 'markvariogram') {
        tmp %<>% mutate(p_val = -1, p_val_adj = -1)
    } else if (de_method == 'smash') {
        tmp %<>% mutate(p_val = -1, p_val_adj = -1, stat = shap_value,
                        cell_type = label)
    } else if (de_method == 'spacgn') {
        tmp %<>% mutate(gene = genes)
    } else if (de_method == 'spatialDE') {
        tmp %<>% mutate(stat = 'LLR', p_val_adj = 0)
    } else if (de_method == 'spatialDE2'){
        tmp %<>% mutate(stat = 'LLR', p_val_adj=0)
    } else if (de_method == 'spagft'){
        tmp %<>% mutate(stat = 'gft_score', p_val_adj=0)
    } else if (de_method == 'spanve'){
        tmp %<>% mutate(stat = 'ent', p_val_adj=0)
    }
    
    if ('type' %in% colnames(tmp)) {
        tmp %<>%
            mutate(de_method = paste0(de_method, '_', type))
    } else {
        tmp %<>%
            mutate(de_method = de_method)
    }
    if (!'stat' %in% colnames(tmp)) tmp %<>% mutate(stat = 0)
    if (!'cell_type' %in% colnames(tmp)) tmp %<>% mutate(cell_type = 'CellType0')
    
    # moransi have p_vals < 0
    tmp %<>% 
        filter(!is.na(p_val))
    
    for (curr_de_method in unique(tmp$de_method)) {
        print(curr_de_method)
        first_de_plot = TRUE
        tmp2 = tmp %>% 
            filter(de_method == curr_de_method) %>%
            dplyr::select(gene, stat, p_val, p_val_adj)
        tmp2$p_val_adj = p.adjust(tmp2$p_val, method = 'BH')
        tmp2$de_binary = tmp2$p_val_adj < 0.05
        
        dat0 = truth %>% 
            dplyr::select(gene, truth) %>% 
            left_join(tmp2 ,by = 'gene'
            ) %>%
            filter(!is.na(truth))
        
        plot_types = c(
            'True positive',
            'True negative',
            'False positive',
            'False negative'
        )
        for (plot_type in plot_types) {
            if (plot_type == 'True positive') {
                genes = dat0 %>%
                    filter(truth == 1, de_binary)
            } else if (plot_type == 'True negative') {
                genes = dat0 %>%
                    filter(truth == 0, !de_binary) 
            } else if (plot_type == 'False positive') {
                genes = dat0 %>%
                    filter(truth == 0, de_binary)
            } else if (plot_type == 'False negative') {
                genes = dat0 %>%
                    filter(truth == 1, !de_binary)
            }
            genes = genes[sample(nrow(genes)),]
            genes_to_plot = genes %>% pull(gene) %>% head(1)
            # empty plot if no genes of this class generated
            if (length(genes_to_plot) < 1) 
                genes_to_plot %<>% c(NA)
            
            meta0 = meta %>%
                dplyr::select(barcode, x, y, label)
            for (gene in genes_to_plot) {
                if (!is.na(gene)) {
                    meta0$expr = expr[gene, meta0$barcode]
                    plot_df = data.frame()
                    conditions = c('Perturbation_1', 'Perturbation_2')
                    for (condition in conditions) {
                        tmp_plot_df = meta0 %>% filter(label == !!condition)
                        fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.015)
                        tmp_plot_df$expr_fit = predict(fit, tmp_plot_df)
                        plot_df %<>% rbind(tmp_plot_df)
                    }
                    set_alpha = 1 # crucial for empty plots
                } else {
                    plot_df = data.frame()
                    for (condition in conditions) {
                        tmp_plot_df = meta0 %>% filter(label == !!condition)
                        tmp_plot_df$expr = NA
                        tmp_plot_df$expr_fit = runif(nrow(tmp_plot_df))
                        plot_df %<>% rbind(tmp_plot_df)
                    }
                    set_alpha = 0 # crucial for empty plots
                }
                
                full_range = range(plot_df$expr_fit)
                brks = c(full_range[1] + 0.1 * diff(full_range), full_range[2] - 0.1 * diff(full_range))
                # full_range = c(NA, quantile(plot_df$expr_fit, probs = 0.99, na.rm = TRUE))
                labels = c('min', 'max')
                
                # color_pal = colorRampPalette(c("yellow", "purple"))(10)
                # color_pal = nr_heat_blue_no_white %>% tail(-5)
                # color_pal = nr_heat_red_no_white %>% tail(-5)
                # color_pal = nr_heat_blue_spatial
                color_pal = pals::kovesi.linear_blue_95_50_c20(100)
                # color_pal = cet_pal(100, name = 'l19') %>% rev()
                
                # plot_df$label = gsub('Perturbation_', 'Label ', plot_df$label)
                plot_df$label = ifelse(plot_df$label == 'Perturbation_1', 'Control', 'Treatment')
                plot_df$de_method = de_method_names$x_name[de_method_names$de_method_name == curr_de_method]
                
                expr_plot = plot_df %>%
                    arrange(-expr_fit) %>%
                    mutate(expr_fit = winsorize(expr_fit, limits = full_range)) %>%
                    ggplot(aes(x = y, y = x, fill = expr_fit)) +
                    rasterise(geom_point(size = 0.1, shape = 21, stroke = 0, alpha = set_alpha), dpi = 600) +
                    scale_color_gradientn(colours = color_pal, name = 'Expression', breaks = brks, labels=labels) +
                    scale_fill_gradientn(colours = color_pal, name = 'Expression', breaks = brks, labels=labels) +
                    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
                    boxed_theme(size_lg = 5, size_sm = 5) +
                    ggtitle(plot_type) +
                    theme(
                        plot.title = element_text(size=5,margin = margin(0,0,-5,0)),
                        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
                        axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.ticks.length.x = unit(0, 'lines'),
                        axis.ticks.length.y = unit(0, 'lines'),
                        legend.position = 'none',
                        panel.spacing = unit(0.05, "cm")
                    ) +
                    coord_fixed()
                expr_plot
                if (first_de_plot) {
                    expr_plot = expr_plot +
                        facet_grid(rows = vars(de_method) ,
                                   cols = vars(label),
                                   switch = 'y'
                        ) + theme(
                            strip.text.y.left = element_text(angle = 0, hjust = 1)
                        )
                    first_de_plot = F
                } else {
                    expr_plot = expr_plot + facet_grid(~label)
                }
                if (de_method != first(de_methods)) {
                    expr_plot = expr_plot +
                        theme(
                            plot.title = element_blank(),
                            strip.text.x = element_blank()
                        )
                }
                expr_plots[[length(expr_plots)+1]] = expr_plot
            }
        }
    }
}
expr_plots[[length(expr_plots)]] = expr_plots[[length(expr_plots)]] +
    theme(
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.14, 'lines'),
        legend.text = element_text(size=5),
        legend.title = element_text(size=5)
    )
p0 = wrap_plots(expr_plots, ncol = 4)
plot_height = 13
plot_width = 12
ggsave(paste0('fig/EFig1/', input, '.pdf'), p0, width = plot_width, height = plot_height, units='cm')

