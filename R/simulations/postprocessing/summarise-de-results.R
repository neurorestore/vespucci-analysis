setwd('~/git/vespucci-analysis/')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(pROC)
library(nebula)
library(pROC)
library(PRROC)
library(mltools)

convert_filename_to_param = function(filename){
	ori_filename = filename
	filename = gsub('.rds', '', basename(filename))
	filename = gsub('.csv', '', filename)
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
	return(values)
}

calculate_de_stats = function(output) {
	auroc_val = 
		tryCatch({
			roc_res = roc(
				predictor = output$de[!is.na(output$de)],
				response = factor(output$truth[!is.na(output$de)], levels = c(0,1))
			)
			auc(roc_res)[1]
		},
		error = function(e){
			message(e)
			return (-1)
		})

	auprc_vals = 
		tryCatch({
			pr_res = pr.curve(scores.class0 = output$de[!is.na(output$de) & output$truth==1], scores.class1 = output$de[!is.na(output$de) & output$truth==0])
			list(
				'auprc_integral' = pr_res$auc.integral,
				'auprc_davis_goadrich' = pr_res$auc.davis.goadrich
			)
		},
		error = function(e){
			message(e)
			return(
				list(
					'auprc_integral' = -1,
					'auprc_davis_goadrich' = -1
				)
			)
		})

	output %<>%
		mutate(
			de_binary_int = as.integer(de_binary)
		)

	tp_genes = output %>%
		filter(de_binary_int == truth, truth == 1) %>%
		pull(gene)
	tn_genes = output %>%
		filter(de_binary_int == truth, truth == 0) %>%
		pull(gene)
	fp_genes = output %>%
		filter(de_binary_int != truth, de_binary == 1) %>%
		pull(gene)
	fn_genes = output %>%
		filter(de_binary_int != truth, de_binary != 1) %>%
		pull(gene)

	tp_size = length(tp_genes)
	tn_size = length(tn_genes)
	fp_size = length(fp_genes)
	fn_size = length(fn_genes)

	sensitivity = tp_size/(tp_size + fn_size)
	specificity = tn_size/(tn_size + fp_size)
	ppv = tp_size/(tp_size + fp_size)
	npv = tn_size/(tn_size + fn_size)

	acc = (tp_size + tn_size)/(nrow(output))

	mcc = mcc(output$de_binary_int, output$truth)

	stats =	data.frame(
		ngenes = nrow(output),
		tp = tp_size,
		tn = tn_size,
		fp = fp_size,
		fn = fn_size,
		auroc = auroc_val,
		auprc_integral = auprc_vals[['auprc_integral']],
		auprc_davis_goadrich = auprc_vals[['auprc_davis_goadrich']],
		sensitivity = sensitivity,
		specificity = specificity,
		ppv = ppv,
		npv = npv,
		acc = acc,
		mcc = mcc
	)
	return (stats)
}

seurat_dir = 'data/simulations/objects/'
outer_de_res_dir = 'data/simulations/de_results/'
de_res_files = c( 
	list.files(paste0(outer_de_res_dir, 'others'), full.names=T),
	list.files(paste0(outer_de_res_dir, 'others_python'), full.names=T)
)
de_res_files = de_res_files[grepl('.rds', de_res_files) | grepl('.csv', de_res_files)]
de_res_files_df = map_df(de_res_files, convert_filename_to_param) %>%
	as_tibble()

parameter_sets = de_res_files_df %>%
	dplyr::select(input, de_prob, seed) %>%
	distinct()

stats_summary_df = data.frame()

for (i in 1:nrow(parameter_sets)) {
	print(i)
	tmp_de_res_files_df = parameter_sets[i,] %>%
		left_join(de_res_files_df)

	parameter_set_prefix = gsub('-de=.*', '', basename(tmp_de_res_files_df$ori_filename[1]))
	sc_file = paste0(seurat_dir, parameter_set_prefix, '.rds')
	sc = readRDS(sc_file)

	gene_features = sc@assays$originalexp@meta.features
	perturb_1_cols = colnames(gene_features)[grepl('Perturbation_1', colnames(gene_features)) & grepl('_is_selected', colnames(gene_features))]
	perturb_2_cols = colnames(gene_features)[grepl('Perturbation_2', colnames(gene_features)) & grepl('_is_selected', colnames(gene_features))]
	gene_features$de_effect = rowSums(gene_features[,c(perturb_1_cols, perturb_2_cols)]) > 0
	gene_features$truth = as.numeric(gene_features$de_effect)

	all_de_methods_genes = c()

	for (de_method in unique(tmp_de_res_files_df$de)){
		print(de_method)
		de_filename = tmp_de_res_files_df %>%
			filter(de == de_method) %>%
			pull(ori_filename)

		if (endsWith(de_filename, '.csv')) {
			tmp = read.csv(de_filename)
		} else if (endsWith(de_filename, '.rds')) {
			tmp = readRDS(de_filename)
		}

		# all manual fixes
		if (de_method == 'squidpy'){
			tmp %<>% mutate(gene = Gene)
		} else if (de_method == 'spacgn'){
			tmp %<>% mutate(gene = genes)
		} else if (de_method == 'spatialDE'){
			tmp %<>% mutate(stat = 'LLR', p_val_adj=p.adjust(p_val))
		}

		if ('type' %in% colnames(tmp)){
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

		post_hoc_grid = tmp %>%
			dplyr::select(cell_type, de_method) %>%
			distinct()
		for (j in 1:nrow(post_hoc_grid)) {
			tmp_row = post_hoc_grid[j,]
			curr_cell_type = tmp_row$cell_type
			curr_de_method = tmp_row$de_method

			tmp_output = tmp %>%
				filter(
					de_method == curr_de_method,
					cell_type == curr_cell_type
				)

			if (length(all_de_methods_genes) == 0) {
				all_de_methods_genes = tmp_output$gene
			} else {
				intersect_genes = intersect(tmp_output$gene, all_de_methods_genes)
				if (length(intersect_genes) > 0 ) {
					all_de_methods_genes = intersect_genes
				} else {
					message(paste0('Cannot find intersection for ', curr_de_method))
				}
			}

			ori_gene_size = nrow(tmp_output[!is.na(tmp_output$p_val),])

			tmp_output %<>% dplyr::select(cell_type, gene, stat, p_val, p_val_adj, de_method)
				
			# somehow moransi have negative p-values
			tmp_output %<>% filter(p_val >= 0)
			tmp_output$p_val_adj = p.adjust(tmp_output$p_val, method = 'BH')
			min_pval = min(tmp_output$p_val[tmp_output$p_val > 0])

			tmp_output %<>%
				mutate(
					de_binary = p_val_adj < 0.05,
					de = ifelse(p_val > 0, -log10(p_val), -log10(min_pval))
				)
			
			tmp_output %<>%
				right_join(gene_features %>% dplyr::select(gene, de_effect, truth), by='gene') %>%
				mutate(
					cell_type = curr_cell_type,
					de_method = curr_de_method,
					de_binary = ifelse(is.na(de_binary), FALSE, de_binary)
				)

			tmp_output %<>% type_convert()
			# compute stats here
			stats = calculate_de_stats(tmp_output)
			stats[is.na(stats)] = 0

			stats = parameter_sets[i,] %>%
				mutate(
					cell_type = curr_cell_type,
					de_method = curr_de_method,
					ori_gene_size = ori_gene_size
				) %>%
				cbind(., stats)
			stats_summary_df %<>% rbind(stats)
		}
	}
}

saveRDS(stats_summary_df, 'data/simulations/summaries/de-results.rds')