library(tidyverse)
library(magrittr)
library(Seurat)
library(hdf5r)
library(Matrix)
library(ggplot2)
library(imager)
library(RNiftyReg)
library(gridExtra)
library(reshape2)
library(sparseMatrixStats)
library(rearrr)
library(rjson)
library(grid)

setwd('~/git/vespucci-analysis/')
source('R/functions/registration_functions.R')

outer_data_dir = 'data/published_data/raw_data/Calcagno2022/'
data_dirs = list.dirs(outer_data_dir, recursive = F)
data_dirs = data_dirs[!endsWith(data_dirs, 'tar')]

folder_to_params = list(
    "V_d1_1" = c('MI', 'd1', 'rep1'),
    "V_d1_2" = c('MI', 'd1', 'rep2'),
    "V_d1_3" = c('MI', 'd1', 'rep3'),
    "V_d7_1" = c('MI', 'd7', 'rep1'),
    "V_d7_2" = c('MI', 'd7', 'rep2'),
    "V_d7_3" = c('MI', 'd7', 'rep3'),
)

sc_list = list()
for (data_dir in data_dirs){
    sample = basename(data_dir)
    if (sample %in% names(folder_to_params)){
        params = folder_to_params[[sample]]
        label = params[1]
        timepoint = params[2]
        replicate = paste0(params, collapse='_')
        
        h5_file = paste0(data_dir, '/filtered_feature_bc_matrix.h5')
        h5 = H5File$new(h5_file)
        print(replicate)
        
        barcodes = readDataSet(h5[['matrix']][['barcodes']])
        barcodes = paste0(replicate, '-', barcodes)
        
        dat_x = readDataSet(h5[['matrix']][['data']])
        # somehow range starts from 0
        row_idx = readDataSet(h5[['matrix']][['indices']]) + 1
        barcodes_idx_map = readDataSet(h5[['matrix']][['indptr']])
        barcodes_idx = unlist(unname(lapply(1:(length(barcodes_idx_map)-1), function(x){
            rep_len = barcodes_idx_map[x+1] - barcodes_idx_map[x]
            rep(x, rep_len)
        })))
        shape = readDataSet(h5[['matrix']][['shape']])
        mat = sparseMatrix(
            i = row_idx,
            j = barcodes_idx,
            x = dat_x,
            dims = shape)
        # print(dim(mat))
        
        genes = readDataSet(h5[['matrix']][['features']][['name']])
        rownames(mat) = genes
        colnames(mat) = barcodes
        
        json_file = paste0(data_dir, '/scalefactors_json.json')
        json_parameters = fromJSON(file=json_file)
        
        tissue_position_file = paste0(data_dir, '/tissue_positions_list.csv')
        tissue_position = read.csv(tissue_position_file, header = F, col.names = c('barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres')) %>%
            mutate(
                barcode = paste0(replicate, '-', barcode),
                replicate = replicate,
                label = label,
                timepoint = timepoint, 
                x = array_col,
                y = array_row,
                low_col = pxl_col_in_fullres * json_parameters$tissue_lowres_scalef, 
                low_row = pxl_row_in_fullres * json_parameters$tissue_lowres_scalef,
                high_col = pxl_col_in_fullres * json_parameters$tissue_hires_scalef, 
                high_row = pxl_row_in_fullres * json_parameters$tissue_hires_scalef,
                low_col_idx = as.integer(low_col),
                low_row_idx = as.integer(low_row),
                high_col_idx = as.integer(high_col),
                high_row_idx = as.integer(high_col)
            ) %>%
            filter(barcode %in% colnames(mat))
        
        rownames(tissue_position) = tissue_position$barcode
        sum(colnames(mat) %in% tissue_position$barcode)
        mat = mat[,tissue_position$barcode]
        
        sc0 = CreateSeuratObject(mat, meta.data = tissue_position)
        sc_list[[length(sc_list)+1]] = sc0   
    }
}

sc = merge(sc_list[[1]], sc_list[2:length(sc_list)])
sc@meta.data %>% dplyr::select(replicate, label, timepoint) %>% distinct()

x_size = 150
y_size = 100
donut_coords = tidyr::crossing(x=1:x_size, y=1:y_size)
start_x = as.integer(median(1:x_size))
start_y = as.integer(median(1:y_size))
r = 50
circle_thickness = 30
r2 = r - circle_thickness
donut_coords %<>%
    mutate(
        in_tissue = ifelse(
            ((x - start_x) * (x - start_x) + (y - start_y) * (y - start_y)) <= r * r, 1, 0),
        in_tissue = ifelse(
            ((x - start_x) * (x - start_x) + (y - start_y) * (y - start_y)) <= r2 * r2, 0, in_tissue)
    )
donut_coords %<>%
    filter(in_tissue == 1)
ori_meta = meta

replicates = unique(meta$replicate)
new_meta = data.frame()
for (i in 1:length(replicates)){
    print(i)
    replicate = replicates[i]
    
    tmp_meta = meta %>%
        filter(replicate == !!replicate)
    tmp_meta = shift_coords(tmp_meta)
    tmp_meta = centralize_coordinates(tmp_meta, x_size = x_size, y_size = y_size)
    new_meta %<>% rbind(tmp_meta)
}

scale_params = data.frame()
for (replicate in unique(meta$replicate)){
    meta0 = meta %>%
        filter(replicate == !!replicate)
    scale_params %<>% rbind(
        data.frame(
            replicate = replicate,
            mean_x = mean(meta0$x),
            sd_x = sd(meta0$x),
            mean_y = mean(meta0$y),
            sd_y = sd(meta0$y)
        )
    )
}

rotation = list(
    'MI_d1_rep1' = c('degree'=90),
    'MI_d1_rep2' = c('degree'=100),
    'MI_d1_rep3' = c('degree'=90),
    'MI_d7_rep1' = c('degree'=250),
    'MI_d7_rep2' = c('degree'=100)
)

filtered_meta = meta %>%
    filter(!replicate %in% names(rotation)) %>%
    mutate(ori_x = x, ori_y = y)
rotated_meta %<>%
    mutate(
        ori_x = x, 
        ori_y = y,
        x = x_rotated,
        y = y_rotated
    )
rotated_meta = rotated_meta[,colnames(filtered_meta)]
meta = rbind(filtered_meta, rotated_meta)

x_size = 75
y_size = 75
new_meta = data.frame()
replicates = unique(meta$replicate)
for (i in 1:length(replicates)){
    print(i)
    replicate = replicates[i]
    
    tmp_meta = meta %>%
        filter(replicate == !!replicate)
    tmp_meta = shift_coords(tmp_meta)
    tmp_meta = centralize_coordinates(tmp_meta, x_size = x_size, y_size = y_size)
    new_meta %<>% rbind(tmp_meta)
}

new_meta = shift_coords(new_meta)
meta = new_meta

x_size = 150
y_size = 120

# barcode, ref_sample, target_sample, new_x, new_y
full_barcodes_coords = data.frame()
full_grid = tidyr::crossing(
    ref_rep = unique(meta$replicate),
    target_rep = unique(meta$replicate)
) %>%
    filter(ref_rep != target_rep)

i = 1
meta$expr = meta$nCount_RNA
for (i in 1:nrow(full_grid)){
    print(i)
    tmp_row = full_grid[i,]
    
    ref_meta = meta %>%
        filter(replicate == tmp_row$ref_rep)
    
    target_meta = meta %>%
        filter(replicate == tmp_row$target_rep)
    
    ref_img = create_bin_img_from_coords(ref_meta, x_size, y_size, interpolate = T, span=0.001)
    plot(as.cimg(ref_img))
    target_img = create_bin_img_from_coords(target_meta, x_size, y_size, interpolate = T, span=0.001)
    plot(as.cimg(target_img))
    
    H = register(target_img, ref_img)
    reg_target_img = applyTransform(H, target_img)
    plot(as.cimg(reg_target_img))
    
    filler_size = x_size * y_size
    filler_image = matrix(1:filler_size, nrow=x_size, ncol=y_size)
    filler_reg_image = applyTransform(H, filler_image, interpolation = 0, nearest = T) %>% as.array()
    
    ## we know there's no zero since it starts from 1
    filler_image_new_coords = melt(filler_reg_image) %>%
        filter(value != 0) %>%
        set_colnames(c('new_x', 'new_y', 'filler_value'))
    
    ## get the barcodes coordinates (i.e. old coords)
    barcodes_coords = target_meta %>%
        dplyr::select(barcode, x, y) %>%
        dplyr::rename(old_x = x) %>%
        dplyr::rename(old_y = y)
    
    filler_image_coords = melt(filler_image) %>%
        set_colnames(c('old_x', 'old_y', 'filler_value')) %>%
        left_join(filler_image_new_coords, by='filler_value', multiple='all')
    
    final_coords = barcodes_coords %>%
        left_join(filler_image_coords, by=c('old_x', 'old_y'), multiple='all') %>%
        filter(!is.na(new_x), !is.na(new_y))
    
    ## finally add to full barcode coords
    full_barcodes_coords %<>%
        rbind(., final_coords %>%
                  mutate(
                      target_sample = tmp_row$target_rep,
                      reference_sample = tmp_row$ref_rep
                  ) %>%
                  dplyr::select(-filler_value) %>%
                  group_by(target_sample, barcode) %>%
                  summarise(
                      final_x = mean(new_x),
                      final_y = mean(new_y)
                  ) %>%
                  mutate(
                      source_sample = tmp_row$ref_rep
                  )
        )
}

mean_barcodes_coords = full_barcodes_coords %>%
    group_by(target_sample, barcode) %>%
    summarise(
        n_count = n(),
        final_x = mean(final_x),
        final_y = mean(final_y)
    ) %>%
    ungroup() %>%
    dplyr::rename(replicate = target_sample)

meta %<>% left_join(mean_barcodes_coords, by=c('replicate', 'barcode'))
meta %<>%
    mutate(
        x = final_x,
        y = final_y
    ) %>%
    dplyr::select(-final_x, -final_y) %>%
    filter(
        !is.na(x),
        !is.na(y)
    )

rownames(meta) = meta$barcode
sc = sc[, rownames(meta)]
meta = meta[colnames(sc),]
all(colnames(sc) == rownames(meta))
sc@meta.data = meta

saveRDS(sc, 'data/published_data/seurat/Calcagno2022.rds')
