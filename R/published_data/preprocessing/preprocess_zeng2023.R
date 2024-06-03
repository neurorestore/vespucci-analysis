library(tidyverse)
library(magrittr)
library(imager)
library(RNiftyReg)
library(gridExtra)
library(reshape2)
library(countsplit)

setwd('~/git/vespucci-analysis/')
source('R/functions/registration_functions.R')

sc = readRDS('data/published_data/raw_data/Zeng2023/Zeng2023_raw.rds')
all(rownames(sc@meta.data) == colnames(sc))

meta = sc@meta.data %>%
    type_convert() %>%
    filter(time == '13months')

scale_down_factor = 100
meta %<>%
    mutate(
        replicate = label.x,
        cell_type = top_level_cell_type
    ) %>%
    group_by(replicate) %>%
    mutate(
        x = as.integer(X/scale_down_factor),
        y = as.integer(Y/scale_down_factor)
    ) %>%
    filter(!is.na(x), !is.na(y)) %>%
    ungroup()

meta %<>%
    mutate(
        marker = ifelse((grepl('CA', cell_type) | cell_type == 'DG'), 'marker', 'background')
    )


ori_meta = meta

replicates = unique(meta$replicate)
i = 1
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

meta = new_meta
rm(new_meta)

square_coords = tidyr::crossing(x=1:x_size, y=1:y_size) %>%
    mutate(in_tissue=1)
square_coords %>%
    ggplot(aes(x=x, y=y)) +
    geom_point(size=1) +
    boxed_theme()

# register shape first
full_barcodes_coords = data.frame()
for (i in 1:length(replicates)){
    print(i)
    replicate = replicates[i]
    
    target_meta = meta %>%
        filter(replicate == !!replicate)
    
    ref_img = create_bin_img_from_coords(square_coords, x_size, y_size)
    plot(as.cimg(ref_img))
    target_img = create_bin_img_from_coords(target_meta, x_size, y_size, interpolate = T, span=0.001)
    plot(as.cimg(target_img))
    
    H = register_shape(target_img, ref_img)
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
                      target_sample = replicate
                  ) %>%
                  dplyr::select(-filler_value) %>%
                  group_by(target_sample, barcode) %>%
                  summarise(
                      final_x = mean(new_x),
                      final_y = mean(new_y)
                  )
        )
}

mean_barcodes_coords = full_barcodes_coords %>%
    group_by(target_sample, barcode) %>%
    summarise(
        n_count = n(),
        x = mean(final_x),
        y = mean(final_y)
    ) %>%
    ungroup() %>%
    dplyr::rename(replicate = target_sample)

meta %<>%
    dplyr::select(-x, -y) %>%
    left_join(mean_barcodes_coords, by=c('replicate', 'barcode'))


full_barcodes_coords = data.frame()
full_grid = tidyr::crossing(
    ref_rep = unique(meta$replicate),
    target_rep = unique(meta$replicate)
) %>%
    filter(ref_rep != target_rep)

i = 1

for (i in 1:nrow(full_grid)){
    print(i)
    tmp_row = full_grid[i,]
    
    ref_meta = meta %>%
        filter(replicate == tmp_row$ref_rep) 
    
    target_meta = meta %>%
        filter(replicate == tmp_row$target_rep)
    
    ref_img = create_img_from_coords(ref_meta, max_x, max_y)
    target_img = create_img_from_coords(target_meta, max_x, max_y)
    plot(as.cimg(ref_img))
    plot(as.cimg(target_img))
        
    H = register(target_img, ref_img)
    reg_target_img = applyTransform(H, target_img)
    plot(as.cimg(reg_target_img))
    
    filler_size = max_x * max_y
    filler_image = matrix(1:filler_size, nrow=max_x, ncol=max_y)
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

marker_colors = c('marker'='black', 'background'='white')

meta %<>%
    dplyr::rename(
        ori_x = x,
        ori_y = y,
        x = final_x,
        y = final_y
    ) %>%
    filter(
        !is.na(x), !is.na(y)
    )
meta %<>% as.data.frame()
rownames(meta) = meta$barcode
    
sc = sc[,rownames(meta)]
sc@meta.data = meta

meta$dens = get_density(meta$x, meta$y)
density_filter = quantile(meta$dens, 0.02)
meta %<>%
    filter(dens >= density_filter)  %>%
    filter(y <= 200)

sc = sc[,rownames(meta)]
saveRDS(sc, 'data/published_data/seurat/Zeng2023.rds')