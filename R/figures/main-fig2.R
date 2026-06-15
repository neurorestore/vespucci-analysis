setwd('~/git/vespucci-analysis/')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(pROC)
library(PRROC)
library(readxl)
# library(lawstat)
# library(nparcomp)
library(gridExtra)
library(cetcolor)
library(ggrastr)
source('R/theme.R')
source('R/functions/utils.R')

#############################################################################-
# Koupourtidou auc
#############################################################################-

sc = readRDS('data/real_data/seurat/Koupourtidou2024.rds')
meta = sc@meta.data %>%
    mutate(barcode = gsub('-', '_', barcode))

auc_res = readRDS('data/real_data/vespucci/Koupourtidou2024-seed=42-nsub=10.rds')[[1]]$aucs
meta %<>% left_join(auc_res)

# interpolate in 2D
fit = loess(auc ~ x * y, data = meta, span = 0.015)
meta$auc_fit = predict(fit, meta)

range = range(meta$auc_fit)
brks = c(range[1] + 0.1 * diff(range), range[2] - 0.1 * diff(range))
labels = format(range, digits = 2)
labels = c(paste0(labels[1], ' '), paste0(' ', labels[2]))

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
    scale_fill_gradientn(colours = auc_pal, name = 'AUC   ', labels = labels, limits = range, breaks = brks) +
    scale_color_gradientn(colours = auc_pal, name = 'AUC   ', labels = labels,limits = range, breaks = brks) +
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
ggsave(paste0("fig/final/Fig2/koupourtidou_AUC.pdf"), p10, width = 2, height = 2.7, units = "cm", useDingbats = FALSE)

###############################################################################-
## Koupourtidou genes ####
###############################################################################-
sc = readRDS('data/real_data/seurat/Koupourtidou2024.rds')
meta = sc@meta.data %>% mutate(barcode = gsub('-', '_', barcode))

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
    brks = c(range[1] + 0.1 * diff(range), range[2] - 0.1 * diff(range))
    labels = c('min', 'max')
    
    pal = nr_heat_red_no_white %>% tail(-5)
    p = plot_df %>%
        mutate(
            label = ifelse(label == 'Intact', 'Intact', '3 dpi'),
            label = factor(label, levels=c('Intact', '3 dpi'))
            ) %>%
        # arrange(-expr) %>%
        ggplot(aes(x = x, y = y, fill = expr)) +
        rasterise(geom_point(size = 0.3,shape = 21, stroke = 0, alpha = 1), dpi = 600) +
        scale_color_gradientn(colours = pal, name = 'Expression',breaks = brks, labels = labels) +
        scale_fill_gradientn(colours = pal, name = 'Expression', breaks = brks, labels = labels) +
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
ggsave('fig/final/Fig2/koupourtidou_genes.pdf', p11, width = 3.5, height = 3.5, units = 'cm', useDingbats = FALSE)

###############################################################################-
## Koupourtidou GO terms ####
###############################################################################-

# load GO
sc = readRDS('data/real_data/seurat_GO/DE/Koupourtidou2024.rds')[[1]]
go_df = readRDS('data/metadata/go_names.rds')
meta = sc@meta.data %>%
    mutate(barcode = gsub('-', '_', barcode))
mat = GetAssayData(sc, slot='counts')

gos_to_plot = c(
    'complement activation, classical pathway',
    'neutrophil activation involved in immune response'
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
    brks = c(range[1] + 0.1 * diff(range), range[2] - 0.1 * diff(range))
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
        rasterise(geom_point(size = 0.3,shape = 21, stroke = 0, alpha = 1), dpi = 600) +
        scale_color_gradientn(colours = pal, name = 'GO module\nscore',breaks = brks, labels = labels) +
        scale_fill_gradientn(colours = pal, name = 'GO module\nscore', breaks = brks, labels = labels) +
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
    if (idx != 2)
        p = p + theme(legend.position = 'none')
    if (idx != 1)
        p = p + theme(strip.text.y.left = element_blank())
    p
    go_plots[[idx]] = p
}
p12 = wrap_plots(go_plots, nrow=1)
ggsave('fig/final/Fig2/koupourtidou_GO-modules.pdf', p12, width = 5, height = 5, units = 'cm', useDingbats = FALSE)

combined_plots = list(expr_plots[[1]], expr_plots[[2]], go_plots[[1]], go_plots[[2]])
p13 = wrap_plots(combined_plots, nrow=1)
ggsave('fig/final/Fig2/koupourtidou_genes-and-GO-modules.pdf', p13, width = 5, height = 5, units = 'cm', useDingbats = FALSE)

de_name_mapping = data.frame(
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

# plot cortical genes
genes_df = read_excel('data/Koupourtidou2024/Belgard_2011_cortical_layer_genes.xls') %>% set_colnames(c('ensembl_id', paste0('L', c('2_3', '4', '5', '6', '6b'),'_uncalib'), 'none_uncalib', paste0('L', c('2_3', '4', '5', '6', '6b'),'_calib'), 'none_calib'))
ensembl_mapping = read.delim('data/Koupourtidou2024/Mus_musculus.NCBIM37.57.gtf.gz', header=F) %>% 
    mutate(
        ensembl_id = str_extract(V9, 'gene_id [^;]+') %>% str_remove('gene_id ') %>% str_remove_all('[\";]'),
        gene = str_extract(V9, 'gene_name [^;]+') %>% str_remove('gene_name ') %>% str_remove_all('[\";]')
    ) %>%
    select(ensembl_id, gene) %>%
    distinct()
genes_df %<>% left_join(ensembl_mapping)

de_methods_df = rbind(
	data.frame(
		res_dir = 'data/real_data/DE_summaries/others/',
		de_method = c('sparkx', 'cside', 'wilcox', 'moransi', 'nnsvg', 'heartsvg', 'scran', 'mast', 'binSpect_kmeans', 'binSpect_rank', 'haystack', 'meringue', 'dCor', 'rv')
	),
	data.frame(
		res_dir = 'data/real_data/DE_summaries/others_python/',
		de_method = c('spacgn', 'squidpy', 'spatialDE', 'spatialDE2', 'spanve', 'spagft')
	)
)
run_grid = tidyr::crossing(
    de_methods_df,
    dataset = 'Koupourtidou2024',
    seed = 42
) %>% 
	mutate(
		input_filename = paste0('data/real_data/seurat/', dataset, '.rds'),
		de_file_suffix = ifelse(res_dir == 'data/real_data/DE_summaries/others_python/', '.csv', '.rds'),
		de_method_corrected = case_when(
			de_method == 'vespucci' ~ '-seed=42-nsub=10-de=nebula_nbgmm',
			startsWith(de_method, 'squidpy') ~ paste0('-de=', gsub('_.*', '', de_method)),
			T ~ paste0('-de=',de_method)
		),
		de_file = 
			paste0(
				res_dir,
				dataset,
				de_method_corrected,
				de_file_suffix
			)
	) %>% filter(file.exists(de_file)) 
full_de_res = data.frame()
for (i in 1:nrow(run_grid)){
    tmp_row = run_grid[i,]
    res = read_de_res(tmp_row$de_file, tmp_row$de_method) %>% dplyr::select(de_method, gene, stat, p_val, p_val_adj)
    full_de_res %<>% rbind(res)
}
ves_res = readRDS('data/real_data/DE_summaries/vespucci/Koupourtidou2024-seed=42-nsub=10-de=nebula_nbgmm.rds') %>% mutate(de_method = 'vespucci', stat=p_val) %>% dplyr::select(de_method, gene, stat, p_val, p_val_adj)
full_de_res %<>% rbind(ves_res)
# full_de_res %<>% group_by(de_method) %>% mutate(pval_rank = rank(p_val, ties='min')) %>% ungroup()
full_de_res %<>% group_by(de_method) %>% mutate(pval_rank = rank(p_val, ties='min'), stat_rank = rank(-abs(stat), tie='min')) %>% ungroup()
layers_to_run = paste0('L', c('2_3', '4', '5'))
layer_genes = map_df(layers_to_run, function(curr_layer) {
		genes_df %>% 
			mutate(layer_p = get(paste0(curr_layer, '_calib'))) %>% 
			arrange(desc(!!sym(paste0(curr_layer, "_calib")))) %>%
			slice_head(n = 1000) %>%
			dplyr::select(gene, layer_p)
	}) %>% group_by(gene) %>% arrange(-layer_p) %>% slice_head(n=1) %>% ungroup()

intersect_df = map_df(unique(full_de_res$de_method), function(de_method){
	filtered_de_res = full_de_res %>% filter(pval_rank <= 50, de_method == !!de_method)
	if (nrow(filtered_de_res) > 50) {
		set.seed(42)
		filtered_de_res %<>% mutate(stat_rank = rank(-abs(stat), tie='min')) %>% arrange(stat_rank)
		cutoff_rank = filtered_de_res$stat_rank[50]
		top_part = filtered_de_res %>% filter(stat_rank < cutoff_rank)
		remaining_slots = 50 - nrow(top_part)
		tied_part = filtered_de_res %>% filter(stat_rank == cutoff_rank)
		filtered_de_res = bind_rows(top_part, tied_part %>% slice_sample(n = remaining_slots))
	}
	data.frame(
		de_method = de_method,
		de_genes_n = nrow(filtered_de_res),
		layer_n = nrow(layer_genes),
		# intersect = sum(layer_genes$gene %in% filtered_de_res$gene),
		# pct_intersect = mean(layer_genes$gene %in% filtered_de_res$gene),
		intersect = sum(filtered_de_res$gene %in% layer_genes$gene),
		pct_intersect = mean(filtered_de_res$gene %in% layer_genes$gene)
	)
})
intersect_df %>% arrange(pct_intersect)

dat0 = intersect_df %>% left_join(de_name_mapping) %>% mutate(
	pct_intersect = pct_intersect * 100,
	label_text = paste0(round(pct_intersect, 2), '%')
) %>% arrange(pct_intersect)


p2 = dat0 %>%
    ggplot(aes(x = reorder(x_name, pct_intersect), y = pct_intersect)) +
    geom_hline(aes(yintercept = 0), size = 0.4, color = 'grey88') +
    geom_segment(aes(xend = x_name, yend = 0), color = 'grey88') +
    geom_point(data = dat0 %>% filter(label_text != "OOT"), shape = 21, stroke = 0.2, size = 0.9, color = 'black', fill = 'grey80') +
	geom_label(aes(label = label_text, y = pct_intersect), color = ifelse(dat0$label_text == 'OOT', 'grey', 'black'), label.padding = unit(0.35, 'lines'), label.size = NA, fill = NA, size = 1.75, show.legend = FALSE, hjust=0)  +
    scale_y_continuous('% of cortical layer genes found as DE', limits = c(0, 65)) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    coord_flip() +
    # boxed_theme() +
    clean_theme() +
    theme(
        axis.title.y = element_blank(),
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.2, 'lines'),
        legend.key.height = unit(0.2, 'lines'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
    )
p2
ggsave('fig/final/Fig2/TBI-pct-intersect-with-layer-genes.pdf', p2, width = 5, height = 6, units = 'cm', useDingbats = FALSE)

# now plot module score plot
sc = readRDS('data/real_data/seurat/Koupourtidou2024.rds')
meta = sc@meta.data
data_dir = 'data/real_data/DE_summaries/check/data/'
data_files = list.files(data_dir, pattern='summ=module', full.names=T)
data_files_df = map_df(data_files, convert_filename_to_params, prefix=F)
data_files_df %<>% filter(dataset == 'Koupourtidou2024', n == 50)

plot_list = list()
for (i in 1:nrow(data_files_df)) {
	x_name = de_name_mapping %>% filter(de_method == data_files_df$de[i]) %>% pull(x_name)
    if (grepl('Squidpy', x_name)) x_name = gsub('Squidpy ', 'Squidpy\n', x_name)
	dat0 = readRDS(data_files_df$ori_filename[i]) 
	print(unique((dat0$n_genes)))
	dat0 %<>% left_join(meta %>% dplyr::select(replicate, label, barcode, x, y))
	conditions = unique(dat0$label)
	plot_df = data.frame()
	for (condition in conditions) {
		tmp_plot_df = dat0 %>% dplyr::filter(label == !!condition)
		fit = loess(expr ~ x * y, data = tmp_plot_df, span = 0.02, degree = 1)
		tmp_plot_df$expr = predict(fit, tmp_plot_df)
		plot_df %<>% rbind(tmp_plot_df)
	}
	range = range(plot_df$expr)
	brks = c(range[1] + 0.1 * diff(range), range[2] - 0.1 * diff(range))
	labels = c('min', 'max')
	pal = nr_heat_red_no_white %>% tail(-5)
	plot_df %<>% mutate(label = factor(ifelse(label == '3dpi', '3 dpi', label), levels=c('Intact', '3 dpi')))
    if (x_name != 'Vespucci') {
        p = plot_df %>%
            # arrange(-expr) %>%
            ggplot(aes(x = x, y = y, fill = expr)) +
            rasterise(geom_point(size = 0.3, shape = 21, stroke = 0, alpha = 1), dpi = 600) +
            scale_color_gradientn(colours = pal, name = 'Expression', breaks = brks, labels = labels) +
            scale_fill_gradientn(colours = pal, name = 'Expression', breaks = brks, labels = labels) +
            guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
            clean_theme(size_lg = 6, size_sm = 5) +
            ggtitle(x_name) +
            theme(
                strip.background = element_blank(),
                strip.text = element_text(size = 4),
                strip.text.y.left = element_text(angle = 0),
                aspect.ratio = 0.6,
                plot.title = element_text(size = 4, margin = margin(0,0,-5,0)),
                plot.margin = unit(c(0.01, 0.01, 0.01, 0.1), "cm"),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y.left = element_line(color = 'black', linewidth=0.1),
                axis.line.x.bottom = element_line(color = 'black', linewidth=0.1),
                legend.position = 'none'
            ) +
            facet_grid(~label) +
            coord_fixed()
            plot_list[[length(plot_list)+1]] = p
    } else {
        p3_1 = plot_df %>%
            # arrange(-expr) %>%
            ggplot(aes(x = x, y = y, fill = expr)) +
            rasterise(geom_point(size = 0.3, shape = 21, stroke = 0, alpha = 1), dpi = 600) +
            scale_color_gradientn(colours = pal, name = 'Expression', breaks = brks, labels = labels) +
            scale_fill_gradientn(colours = pal, name = 'Expression', breaks = brks, labels = labels) +
            guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
            clean_theme(size_lg = 6, size_sm = 5) +
            ggtitle(x_name) +
            theme(
                strip.background = element_blank(),
                strip.text.y.left = element_text(angle = 0),
                aspect.ratio = 0.6,
                plot.title = element_text(size = 5, margin = margin(0,0,-5,0)),
                plot.margin = unit(c(0.01, 0.01, 0.01, 0.1), "cm"),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y.left = element_line(color = 'black', linewidth=0.1),
                axis.line.x.bottom = element_line(color = 'black', linewidth=0.1),
                legend.position = 'none'
            ) +
            facet_wrap(~label, ncol=1) +
            coord_fixed()
    }
}
p3 = plot_grid(p3_1, wrap_plots(plot_list, nrow=2), ncol=2, rel_widths=c(1, 15))
ggsave('fig/final/Fig2/de-genes-module-genes_v0.pdf', p3, width=18, height=2.4, units='cm')
p3 = plot_grid(p3_1, wrap_plots(plot_list, nrow=3), ncol=2, rel_widths=c(1, 7))
ggsave('fig/final/Fig2/de-genes-module-genes_v1.pdf', p3, width=18, height=4, units='cm')
p3 = plot_grid(p3_1, wrap_plots(plot_list, nrow=4), ncol=2, rel_widths=c(1, 5))
ggsave('fig/final/Fig2/de-genes-module-genes_v2.pdf', p3, width=18, height=6, units='cm')

# now bulk de concordance
aucc_df = readRDS('data/Koupourtidou2024/TBI_against_bulk_aucc_df.rds')
plot_df = aucc_df %>% filter(k==2000)
plot_df %<>% left_join(de_name_mapping)
plot_df %<>% mutate(label_text = paste0(round(aucc, 2)))

p4 = plot_df %>%
    ggplot(aes(x = reorder(x_name, aucc), y = aucc)) +
    geom_hline(aes(yintercept = 0), size = 0.4, color = 'grey88') +
    geom_segment(aes(xend = x_name, yend = 0), color = 'grey88') +
    geom_point(data = plot_df %>% filter(label_text != "OOT"), shape = 21, stroke = 0.2, size = 0.9, color = 'black', fill = 'grey80') +
	geom_label(aes(label = label_text, y = aucc), color = ifelse(dat0$label_text == 'OOT', 'grey', 'black'),label.padding = unit(0.35, 'lines'), label.size = NA, fill = NA, size = 1.75, show.legend = FALSE, hjust=0)  +
    scale_y_continuous('% of cortical layer genes found as DE', limits = c(0, 0.16)) +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
    coord_flip() +
    # boxed_theme() +
    clean_theme() +
    theme(
        axis.title.y = element_blank(),
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.2, 'lines'),
        legend.key.height = unit(0.2, 'lines'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
    )
p4
ggsave('fig/final/Fig2/TBI-against-bulk-aucc.pdf', p4, width = 5, height = 6, units = 'cm', useDingbats = FALSE)

# plot ranking
de_methods_df = rbind(
	data.frame(
		res_dir = 'data/real_data/DE_summaries/others/',
		de_method = c('sparkx', 'cside', 'wilcox', 'moransi', 'nnsvg', 'heartsvg', 'scran', 'mast', 'binSpect_kmeans', 'binSpect_rank', 'haystack', 'meringue', 'dCor', 'rv')
	),
	data.frame(
		res_dir = 'data/real_data/DE_summaries/others_python/',
		de_method = c('spacgn', 'squidpy', 'spatialDE', 'spatialDE2', 'spanve', 'spagft')
	)
)
run_grid = tidyr::crossing(
    de_methods_df,
    dataset = 'Koupourtidou2024',
    seed = 42
) %>% 
	mutate(
		input_filename = paste0('data/real_data/seurat/', dataset, '.rds'),
		de_file_suffix = ifelse(res_dir == 'data/real_data/DE_summaries/others_python/', '.csv', '.rds'),
		de_method_corrected = case_when(
			de_method == 'vespucci' ~ '-seed=42-nsub=10-de=nebula_nbgmm',
			startsWith(de_method, 'squidpy') ~ paste0('-de=', gsub('_.*', '', de_method)),
			T ~ paste0('-de=',de_method)
		),
		de_file = 
			paste0(
				res_dir,
				dataset,
				de_method_corrected,
				de_file_suffix
			)
	) %>% filter(file.exists(de_file)) 
full_de_res = data.frame()
for (i in 1:nrow(run_grid)){
    tmp_row = run_grid[i,]
    res = read_de_res(tmp_row$de_file, tmp_row$de_method) %>% dplyr::select(de_method, gene, stat, p_val, p_val_adj)
    full_de_res %<>% rbind(res)
}
ves_res = readRDS('data/real_data/DE_summaries/vespucci/Koupourtidou2024-seed=42-nsub=10-de=nebula_nbgmm.rds') %>% mutate(de_method = 'vespucci', stat=p_val) %>% dplyr::select(de_method, gene, stat, p_val, p_val_adj)
full_de_res %<>% rbind(ves_res)

full_de_res %<>% filter(gene %in% ves_res$gene)
table(full_de_res$de_method)
full_de_res %<>% group_by(de_method) %>% mutate(rank = rank(p_val_adj)) %>% ungroup()
full_de_res %<>% left_join(de_name_mapping)

ves_ranking_df = full_de_res %>% 
    filter(de_method == 'vespucci') %>%
    dplyr::select(gene, rank, x_name)
ranking_df = ves_ranking_df

other_x_names = full_de_res %>% filter(x_name != 'Vespucci') %>% pull(x_name) %>% unique() %>% sort()
for (x_name in other_x_names) {
    rank_check = full_de_res %>% 
        filter(x_name == !!x_name) %>%
        dplyr::rename(new_rank = rank, new_x_name = x_name) %>% 
        inner_join(ves_ranking_df) %>%
        arrange(rank) %>%
        dplyr::select(gene, new_rank, new_x_name) %>% 
        set_colnames(c('gene', 'rank', 'x_name'))
    ranking_df %<>% rbind(rank_check)
}

ranking_df %<>% 
    group_by(x_name) %>% 
    mutate(
        new_rank = (rank - min(rank))/(max(rank) - min(rank))
    ) %>% 
    ungroup()

pal = nr_heat_red_no_white %>% tail(-5)
x_name_order = c('Vespucci', de_name_mapping %>% filter(x_name != 'Vespucci') %>% pull(x_name) %>% sort())
p5 = ranking_df %>% 
    mutate(x_name = factor(x_name, levels=x_name_order)) %>% 
    ggplot(aes(x=x_name, y=1, fill=new_rank)) +
    geom_bar(stat='identity') + 
    scale_fill_gradientn(
        colours = pal, 
        breaks = range(ranking_df$new_rank)) +
    clean_theme() + 
    ggtitle('Gene rank') +
    scale_y_continuous(expand=c(0,0)) +
    ylab('Gene rank based by p-value') +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black', title='Rank')) +
    theme(
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.14, 'lines'),
        legend.text = element_text(size = 5),
        legend.position = 'bottom',
        legend.justification = 'right'
        # legend.title = element_blank()
    )
p5
ggsave('fig/final/Fig2/TBI-gene-rank.pdf', p5, width = 8, height = 8, units = 'cm', useDingbats = FALSE)

p5_0 = ranking_df %>% 
    group_by(x_name) %>% 
    arrange(rank) %>% 
    slice_head(n=1) %>%
    ungroup() %>%
    mutate(x_name = factor(x_name, levels=x_name_order)) %>% 
    ggplot(aes(x=x_name, y=1)) +
    geom_bar(stat='identity') + 
    # scale_fill_gradientn(
    #     colours = pal, 
    #     breaks = range(ranking_df$rank)) +
    clean_theme() + 
    ggtitle('Gene rank') +
    scale_y_continuous(expand=c(0,0)) +
    ylab('Gene rank based by p-value') +
    guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black', title='Rank')) +
    theme(
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.14, 'lines'),
        legend.text = element_text(size = 5),
        legend.position = 'bottom',
        legend.justification = 'right'
        # legend.title = element_blank()
    )
p5_0
ggsave('fig/final/Fig2/TBI-gene-rank-empty.pdf', p5_0, width = 8, height = 8, units = 'cm', useDingbats = FALSE)

