# Generate simulated objects
setwd("~/git/vespucci")
options(stringsAsFactors = FALSE)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-run-rctd.R')
grid = read.delim("sh/grids/published_data/run-rctd.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))

args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(spacexr)
library(Matrix)
library(Seurat)

sp = readRDS(args$sp_input_file)
sc = readRDS(args$sp_input_file)


coords = sp@meta.data %>%
	dplyr::select(x, y)

sp_counts = GetAssayData(sp, slot = 'counts')
nUMI = colSums(sp_counts)
names(nUMI) = rownames(coords)

puck = SpatialRNA(coords, sp_counts, nUMI)

## Prepare single cell ref
sc = readRDS(args$sc_input_file)

meta = sc@meta.data %>%
	dplyr::select(nCount_RNA, nFeature_RNA, cell_type, barcode) %>%
	mutate(
		cell_type = gsub('\\/', '_', cell_type),
		cell_type = gsub(' ', '_', cell_type)
	) %>%
	group_by(cell_type) %>%
	dplyr::filter(n() > 25) %>%
	ungroup()

sc_counts = GetAssayData(sc, slot = 'counts') %>%
  	extract(, meta$barcode)
  # as.matrix()

cell_type_dict = meta %>%
	distinct(cell_type) %>%
	drop_na() %>%
	mutate(cluster_idx = row_number())

cell_type = meta %>%
	left_join(cell_type_dict) %>%
	mutate(cell_type = factor(cluster_idx, levels = cluster_idx, labels = cell_type)) %>%
	pull(cell_type)
names(cell_type) = meta$barcode

nUMI = colSums(sc_counts)
names(nUMI) = meta$barcode

reference = Reference(sc_counts, cell_type, nUMI)

# create RCTD
myRCTD = create.RCTD(puck, reference, 
					max_cores = 10,
					gene_cutoff = args$gene_cutoff,
					gene_cutoff_reg = args$gene_cutoff_reg,
					fc_cutoff = args$fc_cutoff,
					fc_cutoff_reg = args$fc_cutoff_reg,
					UMI_min = args$UMI_min,
					UMI_max = args$UMI_max) ## ran on smallworker

myRCTD = run.RCTD(myRCTD, doublet_mode = "full") ## doublets off for Visium

saveRDS(myRCTD, args$output_filename)