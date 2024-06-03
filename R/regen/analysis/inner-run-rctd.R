# Generate simulated objects
setwd("~/git/vespucci-analysis")
options(stringsAsFactors = FALSE)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-run-rctd.R')
grid = read.delim("sh/grids/regen/analysis/run-rctd.txt")
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

label_check = ifelse(args$label == '7d', 'young', args$label)
if (label_check == 'young' & !'young' %in% sp@meta.data) {
	label_check = 'treated'
}
sp = sp[,sp@meta.data %>% filter(label == label_check) %>% pull(barcode)]
coords = sp@meta.data %>%
	dplyr::select(x, y)

sp_counts = GetAssayData(sp, slot = 'counts')
nUMI = colSums(sp_counts)
names(nUMI) = rownames(coords)

puck = SpatialRNA(coords, sp_counts, nUMI)

## Prepare single cell ref
sc = readRDS(args$sc_input_file)

meta = sc@meta.data %>%
	rename(layer = args$layer) %>%
	dplyr::select(nCount_RNA, nFeature_RNA, layer, barcode) %>%
	group_by(layer) %>%
	dplyr::filter(n() > 25) %>%
	ungroup()

sc_counts = GetAssayData(sc, slot = 'counts') %>%
  	extract(, meta$barcode)
  # as.matrix()

cell_type_dict = meta %>%
	distinct(layer) %>%
	drop_na() %>%
	mutate(cluster_idx = row_number())

cell_type = meta %>%
	left_join(cell_type_dict) %>%
	mutate(layer = factor(cluster_idx, levels = cluster_idx, labels = layer)) %>%
	pull(layer)
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