setwd('~/git/vespucci-analysis')
library(tidyverse)
library(magrittr)
library(Matrix)
library(nebula)
library(pROC)
library(PRROC)
library(mltools)
library(fastglm)

get_spatial_acc = function(meta0, auc_vals, method){
	meta = meta0 %>% 
		dplyr::select(barcode, label_sim, de_bool) %>%
		left_join(
			auc_vals,
			by = 'barcode'
		)	
	meta %<>%
		filter(!is.na(auc))
	auroc_val = 
		tryCatch({
			roc_res = roc(
				predictor = meta$auc,
				response = factor(meta$de_bool, levels = c(0,1))
			)
			auc(roc_res)[1]
		},
		error = function(e){
			message(e)
			return (-1)
		})

	auprc_vals = 
		tryCatch({
			pr_res = pr.curve(scores.class0 = meta$auc[meta$de_bool==1], scores.class1 = meta$auc[meta$de_bool==0])
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
	out_df = tmp_row %>%
		dplyr::select(-ori_filename, -ves_filename, -mag_filename) %>%
		tidyr::crossing(
			data.frame(
				'method' = method,
				'val' = c(auroc_val, auprc_vals$auprc_integral),
				'type' = c('AUROC', 'AUPRC')
			)
		)
	return(out_df)
}

data_files = list.files('data/simulations/objects/', full.names=T)
ves_dir = 'data/simulations/vespucci/'
ves_dir = 'data/simulations/brute_force/'

out_grid = data.frame(
		input_file = data_file
	) %>%
	mutate(
		obj = basename(input_filename),
		input = gsub('de-prob.*='),
		vespucci_file = paste0(ves_dir, obj),
		brute_force_file = paste0(brute_force_dir, obj)
	)

output_df = data.frame()

for (i in 1:nrow(out_grid)) {
	tmp_row = out_grid[i,]
	print(tmp_row)
	input_file = tmp_row$input_file

	sc = readRDS(input_file)	
	meta0 = sc@meta.data %>%
		mutate(
			barcode = Cell
		)
	
	labels_to_run = unique(meta0$label_sim)
	labels_to_run = labels_to_run[labels_to_run != 'label1']

	for (label in labels_to_run) {
		if (label == 'background') {
			meta0$de_bool = as.numeric(meta0$label_sim != 'background')
			truth = 'non_background'
		} else {
			meta0$de_bool = as.numeric(meta0$label_sim == label)
			truth = label
		}

		brute_force_res = readRDS(tmp_row$brute_force_file) %>%
			type_convert() %>%
			dplyr::select(barcode, auc)
		out_df = get_spatial_acc(meta0, mag_auc, 'brute_force') %>%
			mutate(
				truth = truth
			)
		output_df %<>% rbind(out_df)

		ves_auc = readRDS(tmp_row$vespucci_file)$spatial_auc_result$aucs
		out_df = get_spatial_acc(meta0, ves_auc, 'vespucci') %>%
			mutate(
				truth = truth
			)
		out_df = out_df[,colnames(output_df)]
		output_df %<>% rbind(out_df) 
	}
}

saveRDS(output_df, 'data/simulations/summaries/spatial-acc-summary.rds')
