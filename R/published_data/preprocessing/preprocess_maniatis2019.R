library(tidyverse)
library(magrittr)
library(Matrix)
library(Seurat)

setwd('~/git/vespucci-analysis/')
source('R/functions/registration_functions.R')

mapping_metadata = read.csv('metadata/Maniatis_GSE120374_raw_to_reg_metadata.csv', row.names=1)
mapping_metadata %<>%
	mutate(
		raw_barcode = paste0(raw_x, '_', raw_y)
	)

data_dir = 'data/published_data/Maniatis_2019/raw_count_data/'
data_files = list.files(data_dir, full.names=T)
data_files = data_files[!grepl('\\.jpg\\.gz', data_files)]

annotation_files = data_files[grepl('annotations', data_files)]
counts_files = data_files[grepl('aligned_counts_ID', data_files)]

full_metadata = data.frame()
full_count_df = data.frame()

for (annotation_file in annotation_files){
	id_name = gsub('\\.tsv.*', '' , basename(annotation_file))
	
	annotation_df = read.table(annotation_file, header=T)
	
	annotation_df %<>%
		mutate(
			raw_barcode = paste0(xPos, '_', yPos)
		)
	rownames(annotation_df) = annotation_df$raw_barcode
	
	count_files_with_id = counts_files[grepl(id_name, counts_files)]

	for (count_file in count_files_with_id){
		GSM_code = gsub('_stdata.*', '', basename(count_file))
		count_data = read.table(count_file, row.names=NULL, check.names=F)
		genes = as.character(count_data[,1])
		count_data = Matrix(as.matrix(count_data[,-1]), sparse=T)

		rownames(count_data) = genes

		duplicated_genes = unique(genes[duplicated(genes)])
		duplicated_genes_idxs = which(genes %in% duplicated_genes)

		duplicated_genes_counts = do.call(rbind, map(duplicated_genes, function(x){
			duplicated_gene_idxs = which(genes == x)
			duplicated_genes_idxs = c(duplicated_genes_idxs, duplicated_gene_idxs)
			gene_count_data = count_data[duplicated_gene_idxs,]
			gene_count_data = colSums(gene_count_data)
		}))
		rownames(duplicated_genes_counts) = duplicated_genes

		count_data = count_data[-duplicated_genes_idxs,]
		count_data = rbind(count_data, duplicated_genes_counts)
		count_data = count_data[,colnames(count_data) %in% annotation_df$raw_barcode]

		temp_metadata = annotation_df[colnames(count_data),] %>%
			mutate(
				GSM_code = GSM_code,
				barcode = paste0(GSM_code, '_', raw_barcode)
			)

		temp_metadata %<>%
			left_join(mapping_metadata, by=c('GSM_code', 'raw_barcode'))
		rownames(temp_metadata) = temp_metadata$barcode
		full_metadata = rbind(full_metadata, temp_metadata)

		stopifnot(all(colnames(count_data) == temp_metadata$raw_barcode))
		colnames(count_data) = temp_metadata$barcode

		summ = summary(count_data)
		count_df = data.frame(
			gene = rownames(count_data)[summ$i],
			barcode = colnames(count_data)[summ$j],
			count = summ$x
		)
		full_count_df = rbind(full_count_df, count_df)
	}
}

combined_mat = xtabs(count~gene+barcode, full_count_df, sparse=T)

sum(!colnames(combined_mat) %in% rownames(full_metadata))
sum(!rownames(full_metadata) %in% colnames(combined_mat))

dim(combined_mat)
full_metadata = full_metadata[colnames(combined_mat),]

sum(is.na(full_metadata$registered_x))
sc = CreateSeuratObject(counts = combined_mat, meta.data = full_metadata) 

## Meta data downloaded from: https://als-st.nygenome.org/
meta = read_delim("data/published_data/raw/maniatis2019/metadata/mouse_sample_names_sra.tsv") %>%
	dplyr::rename(label = breed, timepoint = age)

sc_meta = sc@meta.data %>%
	separate(GSM_code, c("gsm", "x1", "x2", "x3"), "_") %>%
	mutate(replicate = ifelse(!is.na(x3), paste0(x1, "_", x2, "_", x3),
							paste0(x1, "_", x2))) %>%
	mutate(barcode = colnames(sc)) %>%
	dplyr::select(barcode, orig.ident, nCount_RNA, nFeature_RNA, replicate, 
				registered_x, registered_y) %>%
	# flip 90 degrees
	dplyr::rename(x = registered_y, y = registered_x) %>%
	# mirror horizontal
	mutate(x = -x, y = -y) %>%
	# filter some poorly registered barcodes
	filter(y > -6.5, y < 6.5, x > -7, x < 7) %>%
	drop_na() %>%
	left_join(meta %>% dplyr::select(replicate, label, timepoint) %>% distinct()) %>%
	set_rownames(.$barcode)

# subset the data now
expr = GetAssayData(sc, assay = 'RNA', slot = 'counts') %>%
	extract(, sc_meta$barcode)

# create new seurat object
new_sc = CreateSeuratObject(expr, meta.data = sc_meta) %>%
	NormalizeData()

# test some plots
pdat = data.frame(value = new_sc@assays$RNA@data['Aif1', ]) %>%
	bind_cols(sc_meta) %>%
	arrange(value)
  
# remove the outliers
coords = new_sc@meta.data %>% dplyr::select(x, y) %>%
	dplyr::rename(z = x) %>%
	mutate(x = 1)

outliers = edge_detect(x = coords$x, y = coords$y, z = coords$z, 
						barcodes = colnames(new_sc), peakdist = 20,
						span = 0.3, nbin = 4)

pdat$outlier = ifelse(pdat$barcode %in% outliers$outliers, 1, 0)

# define the barcodes to keep
keep = pdat %>% filter(outlier == 0) %>% pull(barcode)

# subset and save this final object
final_sc = subset(new_sc, cells = keep)

# save
saveRDS(final_sc, "data/published_data/seurat/Maniatis2019.rds")