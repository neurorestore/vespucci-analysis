library(tidyverse)
library(magrittr)
library(imager)
library(RNiftyReg)
library(gridExtra)
library(reshape2)

shift_coords = function(tmp_meta) {
    ## make sure it starts from 0
    min_x = min(tmp_meta$x)
    x_size = max(tmp_meta$x)
    min_y = min(tmp_meta$y)
    y_size = max(tmp_meta$y)
    tmp_meta %<>% mutate(x = x-min_x, y = y-min_y)
    return (tmp_meta)
}

centralize_coordinates = function(tmp_meta, x_size, y_size) {
    shift_x = floor(median(1:x_size))
    shift_y = floor(median(1:y_size))
    tmp_meta %<>%
        mutate(
            ori_x = x,
            ori_y = y,
            x = x - mean(ori_x),
            y = y - mean(ori_y),
            x = x + shift_x,
            y = y + shift_y,
            x = as.integer(x),
            y = as.integer(y),
        )
    return (tmp_meta)
}

create_bin_img_from_coords = function(tmp_meta, x_size, y_size, interpolate = F, span=0.03) {
    if (interpolate == T){
        frame = tidyr::crossing(x=1:x_size,y=1:y_size) %>%
            left_join(tmp_meta %>% mutate(value = 255) %>% dplyr::select(x, y, value), by = c('x', 'y')) %>%
            mutate(value = ifelse(is.na(value), 0, 255))
        fit = loess(value ~ x * y, data = frame, span = span, degree = 1)
        output_image = tidyr::crossing(x = 1:x_size, y = 1:y_size) %>% 
            mutate(value = predict(fit, newdata = .)) %>% 
            mutate(value = ifelse(is.na(value), 0, value)) %>%
            mutate(value = round(value)) %>%
            as.cimg()
    } else {
        output_image = matrix(0, nrow=x_size, ncol=y_size)
        for (j in 1:nrow(tmp_meta)){
            output_image[tmp_meta$x[j], tmp_meta$y[j]] = 255
        }   
        output_image = as.cimg(output_image)
    }
    # plot(output_image)
    return(output_image)
}

register_shape = function(source, target) {
    # Compute non linear deformation: "warping"
    H <- forward(niftyreg.nonlinear(source, target,
                                    symmetric = F, nBins = 16, maxIterations = 800,
                                    bendingEnergyWeight = 5e-5, nLevels = 10), verbose = T)
    return(H)
}


create_img_from_coords = function(tmp_meta, max_x, max_y) {
    output_image = matrix(255, nrow=max_x, ncol=max_y)
    marker_meta = tmp_meta %>%
        filter(marker != 'background')
    for (j in 1:nrow(marker_meta)){
        output_image[marker_meta$x[j], marker_meta$y[j]] = 0
    }
    
    # output_image = as.cimg(output_image)
    return(output_image)
}

register = function(source, target) {
    # Initialise registeration: linear deformation
    init <- forward(niftyreg.linear(source, target))
    # Compute non linear deformation: "warping"
    H <- forward(niftyreg.nonlinear(source, target, init=init,
                                    symmetric = F, nBins = 16, maxIterations = 800,
                                    bendingEnergyWeight = 5e-5, nLevels = 10), verbose = T)
    return(H)
}