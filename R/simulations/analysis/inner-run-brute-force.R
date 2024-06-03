setwd("~/git/vespucci-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'time-full-magellan.R')
grid = read.delim("sh/grids/simulations/time-full-magellan.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))

args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(Seurat)
library(Augur)
library(Matrix)
library(RANN)
library(Vespucci)

# control seeding
set.seed(args$seed)

# load in the object
sc = readRDS(args$input_file)

# get the input objects
input = GetAssayData(sc, slot = 'counts')
meta = sc@meta.data

nn = get_nn(meta)
input %<>% Augur::select_variance()

barcodes = colnames(input)

output_df = data.frame()
for (barcode in barcodes){

	# grab the barcodes to keep
	barcodes_keep = nn %>%
		dplyr::filter(source == barcode) %>%
		group_by(target_group) %>%
		dplyr::filter(rank %in% seq(1,args$k)) %>%
		pull(target)

	# subset the object and prepare for Augur
	input0 = input %>% extract(, barcodes_keep)
	meta0 = meta %>% extract(barcodes_keep, ) %>% set_rownames(barcodes_keep)

	meta0$cell_type = 'brute-force'

	# run Augur
	augur_barcode_start_time = Sys.time()
	print('Running Augur... ')
	augur = calculate_auc(input = input0,
						meta = meta0,
						n_subsamples = args$n_subsamples,
						var_quantile = 1,
						n_threads = 1)
	print('Done')
	augur_barcode_end_time = Sys.time()

	output_df %<>% rbind(., data.frame(
		'barcode'=barcode,
		'auc'=augur$AUC$auc,
		'augur_barcode_start_time' = augur_barcode_start_time,
		'augur_barcode_end_time' = augur_barcode_end_time,
		'augur_barcode_time' = augur_barcode_end_time - augur_barcode_start_time
	))
}

saveRDS(output_df, args$output_file)