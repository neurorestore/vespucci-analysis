setwd('~/git/vespucci-analysis/')
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
library(countsplit)

source('R/theme.R')
source('R/functions/registration_functions.R')

data_dir = 'data/published_data/raw_data/Koupourtidou2024/'
samples = c('GSM7068162_D', 'GSM7068163_A')

sc_list = list()
for (sample in samples) {
    label = ifelse(sample == 'GSM7068162_D', 'Intact', '3dpi')
    h5_file = paste0(data_dir, sample, '_filtered_feature_bc_matrix.h5')
    
    mat = Read10X_h5(h5_file)
    colnames(mat) = paste0(label, '-', colnames(mat))
    
    json_file = paste0(data_dir, sample, '_scalefactors_json.json')
    json_parameters = fromJSON(file=json_file)
    
    tissue_position_file = paste0(data_dir, sample, '_tissue_positions_list.csv')
    tissue_position = read.csv(tissue_position_file, header = F, col.names = c('barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres')) %>%
        mutate(
            barcode = paste0(label, '-', barcode),
            # replicate = replicate,
            label = label,
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

sc = merge(sc_list[[1]], sc_list[[2]])
meta = sc@meta.data
meta %>%
    ggplot(aes(x=x, y=y)) +
    geom_point(size=0.1) +
    facet_wrap(~label) +
    boxed_theme()

meta1 = meta %>%
    filter(label == '3dpi')
split_line = data.frame(
    x = meta1$x
) %>%
    mutate(
        y_line = 43 + (-0.1*x)
    ) %>%
    distinct()
new_meta1 = data.frame()
x_to_test = sort(unique(meta1$x))
for (curr_x in x_to_test) {
    print(curr_x)
    tmp_meta1 = meta1 %>%
        filter(x == curr_x) %>%
        mutate(replicate= label)
    print(paste0('# of rows: ', nrow(tmp_meta1)))
    y_val = split_line %>% filter(x == curr_x) %>% pull(y_line)
    for (j in 1:nrow(tmp_meta1)) {
        tmp_row = tmp_meta1[j,]
        if (tmp_meta1$y[j] > y_val) {
            tmp_meta1$replicate[j] = paste0(tmp_meta1$replicate[j], '_rep1')
        } else {
            tmp_meta1$replicate[j] = paste0(tmp_meta1$replicate[j], '_rep2')
        }
    }
    new_meta1 %<>% rbind(tmp_meta1)
}


new_meta = rbind(new_meta1, new_meta2)
new_meta = new_meta[colnames(sc),]
sc@meta.data = new_meta

meta = sc@meta.data
meta$replicate = factor(meta$replicate, levels=c('Intact_rep1', 'Intact_rep2', '3dpi_rep1', '3dpi_rep2'))

expr0 = GetAssayData(sc, slot='counts') %>% NormalizeData()
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
    'Intact_rep1' = c('degree'=168),
    '3dpi_rep2' = c('degree'=190)
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

flipped = c(
    'Intact_rep1', '3dpi_rep1'
)

flipped_meta = data.frame()
meta$expr = expr0[gene, rownames(meta)]
for (replicate in flipped){
    meta0 = meta %>% filter(replicate == !!replicate) %>% data.frame()
    scale_param = scale_params %>% filter(replicate == !!replicate)
    
    meta0 %<>% 
        mutate(
            scaled_x = (x - scale_param$mean_x)/scale_param$sd_x,
            scaled_y = (y - scale_param$mean_y)/scale_param$sd_y,
            x_flipped_scaled = -scaled_x,
            y_flipped_scaled = scaled_y,
            x_flipped = as.integer((x_flipped_scaled * scale_param$sd_x) + scale_param$mean_x),
            y_flipped = as.integer((y_flipped_scaled * scale_param$sd_y) + scale_param$mean_y)
        )
    flipped_meta %<>% rbind(meta0)   
}

filtered_meta = meta %>%
    filter(!replicate %in% flipped)
flipped_meta %<>%
    mutate(
        x = x_flipped,
        y = y_flipped
    )
flipped_meta = flipped_meta[,colnames(filtered_meta)]
meta = rbind(filtered_meta, flipped_meta)

new_meta = data.frame()
for (replicate in unique(meta$replicate)){
    meta0 = meta %>%
        filter(replicate == !!replicate)
    meta0$x = scale(meta0$x)[,1]
    meta0$y = scale(meta0$y)[,1]
    new_meta %<>% rbind(meta0)
}

# filter accordingly
meta = new_meta %>%
    mutate(expr = nCount_RNA) %>%
    filter(
        !(replicate == 'Intact_rep1' & x < -1),
        !(replicate == 'Intact_rep1' & y < -1.4),
        !(replicate == 'Intact_rep2' & x < -1),
        !(replicate == 'Intact_rep2' & y < -1.4),
        !(replicate == '3dpi_rep1' & y < -1.5),
        !(replicate == '3dpi_rep1' & x < -1.4),
        !(replicate == '3dpi_rep1' & x > 2),
        !(replicate == '3dpi_rep2' & x < -1.2),
    )

replicates = unique(meta$replicate)

meta$x = meta$x * 100
meta$y = meta$y * 100

x_size = 200
y_size = 200
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


new_meta = shift_coords(new_meta)
ori_meta = meta
meta = new_meta

x_size = max(meta$x)
y_size = max(meta$y)

# barcode, ref_sample, target_sample, new_x, new_y
full_barcodes_coords = data.frame()
full_grid = tidyr::crossing(
    ref_rep = unique(meta$replicate),
    target_rep = unique(meta$replicate)
) %>%
    filter(ref_rep != target_rep)
meta$expr = meta$nCount_RNA

for (i in 1:nrow(full_grid)){
    print(i)
    tmp_row = full_grid[i,]
    
    ref_meta = meta %>%
        filter(replicate == tmp_row$ref_rep)
    
    target_meta = meta %>%
        filter(replicate == tmp_row$target_rep)
    
    ref_img = create_img_from_coords(ref_meta, x_size, y_size)
    target_img = create_img_from_coords(target_meta, x_size, y_size)
    plot(as.cimg(ref_img))
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
                  ungroup() %>%
                  mutate(
                      source_sample = tmp_row$ref_rep,
                      section = tmp_row$target_section
                  )
        )
}

final_meta = full_barcodes_coords %>%
    group_by(target_sample, barcode) %>%
    summarise(mean_x = mean(final_x), mean_y = mean(final_y)) %>%
    ungroup() %>%
    dplyr::rename(replicate = target_sample)
final_meta %<>% left_join(meta %>% dplyr::select(barcode, expr))

meta = ori_meta %>%
    dplyr::select(-x, -y) %>%
    left_join(
        final_meta %>%
            dplyr::select(barcode, mean_x, mean_y) %>%
            dplyr::rename(x = mean_x, y = mean_y)
    ) %>%
    filter(!is.na(x), !is.na(y))

x_rad = 350
y_rad = 350
angles = seq(0, 2*pi, 0.01) 
arc_coords = 
    data.frame(
        angle = angles,
        x = x_rad * cos(angles),
        y = y_rad * sin(angles)
    ) %>%
    filter(x > 0, y > 0) %>%
    arrange(x)
arc_coords %>% ggplot(aes(x=x,y=y)) + geom_point()

filtered_meta = data.frame()
for (i in 2:nrow(arc_coords)) {
    arc_coord = arc_coords[i,]
    prev_arc_coord = arc_coords[i-1,]
    tmp_meta = meta %>%
        filter(x < arc_coord$x, x > prev_arc_coord$x, y < arc_coord$y)
    # if (nrow(tmp_meta) > 1) break
    filtered_meta %<>% rbind(tmp_meta)
}

rownames(filtered_meta) = filtered_meta$barcode

meta = filtered_meta
meta$replicate = as.character(meta$replicate)
mat = GetAssayData(sc, slot='count')
mat = mat[,rownames(meta)]

sc = CreateSeuratObject(mat, meta.data=meta)
saveRDS(sc, 'data/published_data/seurat/Koupourtidou2024.rds')
