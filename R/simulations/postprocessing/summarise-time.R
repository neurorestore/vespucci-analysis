setwd("~/git/vespucci")

library(tidyverse)
library(magrittr)

brute_force_dir = 'data/simulations/brute_force/'
vespucci_files = list.files('data/simulations/vespucci_results/', full.names=T)

time_df = data.frame()
for (vespucci_file in vespucci_files){
	obj = basename(vespucci_file)
	obj = gsub('\\.rds', '', obj)

	ves_res = readRDS(vespucci_file)
	time_tracking = ves_res$time_tracking
	vespucci_time = 
		time_tracking$global$time +
		sum(time_tracking$auc$time) +
		sum(time_tracking$cor_convergence$model_time) +
		sum(time_tracking$cor_convergence$predict_time)

	brute_force_file = paste0(brute_force_dir, obj, '.rds')
	brute_force_time = sum(readRDS(brute_force_file)$augur_barcode_time)

	input = gsub('-.*', '', obj)
	input = gsub('input=', '', input)

	time_df %<>% rbind(
		'obj' = obj,
		'input' = input,
		'time' = vespucci_time,
		'method' = 'vespucci'
	)
	time_df %<>% rbind(
		'obj' = obj,
		'input' = input,
		'time' = brute_force_time,
		'method' = 'brute_force'
	)
}

saveRDS(time_df, 'data/simulations/summaries/time-summary.rds')