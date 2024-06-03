setwd("~/git/vespucci-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-generate-simulations.R')
grid = read.delim("sh/grids/simulations/test-new-splatter.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),type = typeof(grid[[param_name]]))

args = parser$parse_args()
print(args)

library(nebula) ## somehow only version 1.4.1 works
library(pROC)
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
source("R/functions/simulation_functions.R")
source("R/functions/splatter_modified_functions.R")
library(scater)
library(countsplit)

set.seed(args$seed)

params = load_spatial_parameters()	
false_DE = FALSE

if (args$input_type == 'circle') {
	multi_de_list = list(
		'Perturbation' = list(
			de_prob = args$de_prob,
			de_effect_size = c(0, args$de_effect_size),
			de_fac_scale = 0.4,
			de_down_prob = 0.5,
			de_selected_idx = c(2),
			de_direction = 'both',
			rep_de_facloc_jitter = c(0.3)
			)
	)
} else if (args$input == 'circle_overlap') {
	multi_de_list = list(
		'Perturbation' = list(
			de_prob = args$de_prob,
			de_effect_size = c(0, 2, 0.5),
			de_fac_scale = 0.4,
			de_down_prob = 0.5,
			de_selected_idx = c(2, 3),
			de_direction = 'both',
			rep_de_facloc_jitter = c(0, 0.1, 0.1)
		)
	)
} else if (args$input == 'flag') {
	multi_de_list = list(
		'Perturbation' = list(
			de_prob = args$de_prob,
			de_effect_size = c(0, 0.4, 0.8, 1.2, 2),
			de_fac_scale = 0.4,
			de_down_prob = 0.5,
			de_selected_idx = c(2, 3, 4, 5),
			de_direction = 'both',
			rep_de_facloc_jitter = c(0, 0.1, 0.1, 0.1, 0.1)
		)
	)
} else if (args$input == 'stripes') {
	multi_de_list = list(
		'Perturbation' = list(
			de_prob = args$de_prob,
			de_effect_size = c(0, 0.5, 2, 0.5),
			de_fac_scale = 0.4,
			de_down_prob = 0.5,
			de_selected_idx = c(2, 3, 4),
			de_direction = 'both',
			rep_de_facloc_jitter = c(0, 0.1, 0.1, 0.1)
		)
	)
} else if (args$input == 'false') {
	multi_de_list = list(
		'Perturbation' = list(
			de_prob = args$de_prob,
			de_effect_size = 2,
			de_fac_scale = 0.4,
			de_down_prob = 0.5,
			de_selected_idx = 1:6,
			de_direction = args$de_direction,
			rep_de_facloc_jitter = c(0)
		)
	)
	false_DE = TRUE
}

sim_grid = simulated_grid(n_cols=100) %>%
	apply_spatial_perturbation(
		perturbation = args$input_type
	)

# only works for binary label for now
sim_grid$Cell = paste0('Cell', 1:10000)

# remove background here
if (args$input != 'false') {
	de_specific_cells_list =
		list(
			'Perturbation' = sim_grid %>%
				filter(label != 'background') %>%
				pull(Cell)
		)

	perturb_labels = vapply(sim_grid$label, function(x) {
		ifelse(x == 'background', 
			'Perturbation_1',
			gsub('label', 'Perturbation_', x)
		)
		},
		as.character(1)
	)
	perturb_factor_levels = paste0('Perturbation_', sort(as.numeric(unique(gsub('.*_', '', perturb_labels)))))
} else {
	de_specific_cells_list = NULL
	perturb_labels = sample(paste0('Perturbation', 1:6), nrow(sim_grid), replace=T)
	perturb_factor_levels = unique(perturb_labels)
}

sim_grid$Perturbation = perturb_labels

# no cell type effect
sim_grid %<>%
	mutate(
		imagecol = jitter(imagecol, amount = 1), 
		imagerow = jitter(imagerow, amount = 2),
		cell_type = 'CellTypeB'
	)
# must match multi_de_list

preset_sim_metadata = sim_grid
preset_sim_metadata$Perturbation = factor(preset_sim_metadata$Perturbation, levels = perturb_factor_levels)
n_reps_per_group = args$n_reps_per_group
rep_depth_jitter = args$rep_depth_jitter

params@batch.rmEffect = FALSE

update_list = list(
	batchCells = batchCells,
	group.prob = group.prob,
	seed = 1
)
params = setParams(params, update=update_list)
params = splatter:::expandParams(params)

# Set random seed
seed = getParam(params, "seed")
set.seed(seed + args$seed)

sc = splat_simulate(
	params = params,
	preset_sim_metadata = preset_sim_metadata,
	multi_de_list = multi_de_list,
	de_specific_cells_list = de_specific_cells_list,
	seed = seed,
	false_DE = false_DE
)

saveRDS(sc, args$object_output_filename)