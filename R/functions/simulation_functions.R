# Setup some functions to clean up the simulation code
# setwd("~/git/vespucci")
library(tidyverse)
library(magrittr)
library(Matrix)

# load parameters derived from our own sc-walk data
load_spatial_parameters = function() {
	params = readRDS("data/simulations/parameters_spatial.rds")
	return(params)
}

# function to setup a 'Spatial' simulated grid
simulated_grid = function(n_rows = 100, n_cols = 100) {
	coords = matrix(1, nrow = n_rows, ncol = n_cols) %>%
		as.data.frame() %>%
		mutate(X = row_number()) %>%
		gather(Y, value, -X) %>%
		mutate(Y = gsub("V", "", Y)) %>%
		type_convert() %>%
		# now, re-scale to put this in the same coordinate space/sizing as the
		# main spatial object in sc-walk. Note, this is hard-coded in.
		mutate(
			Y = scales::rescale(Y, c(25, 336)), 
			X = scales::rescale(X, c(27, 482))) %>%
		dplyr::rename(imagecol = X, imagerow = Y) %>%
		mutate(tissue = 1) %>%
		mutate(row = imagerow, col = imagecol) %>%
		dplyr::select(tissue, row, col, imagerow, imagecol)
}

apply_spatial_perturbation = function(
	sim, perturbation = 'circle',  
	circle_radius = 100, # only for circle
	circle_thickness = 40 # only for circle
	) {
	if (perturbation == 'circle') {
		start_x = 254
		start_y = 180
		r = circle_radius
		r2 = r - circle_thickness
		sim %<>%
			mutate(label = ifelse(
				((imagecol - start_x) * (imagecol - start_x) +
						(imagerow - start_y) * (imagerow - start_y)) <= r * r,
				'perturbation', 'background')) %>%
			mutate(label = ifelse(
				((imagecol - start_x) * (imagecol - start_x) +
						(imagerow - start_y) * (imagerow - start_y)) <= r2 * r2,
				'background', label))		
		labels = c('label1', 'label2')
		labels = vapply(sim$label, function(x) {
			ifelse(x != 'background', sample(labels, 1), x)
		},
			as.character(1)
		)
		sim$label = labels
	} else if (perturbation == 'circle_overlap') {
		start_x1 = 214
		start_x2 = 294
		start_y = 180
		r = 100
		sim %<>%
			mutate(label = ifelse(
				((imagecol - start_x1) * (imagecol - start_x1) +
					(imagerow - start_y) * (imagerow - start_y)) <= r * r,
				'label2', 'background')) %>%
			mutate(label = ifelse(
				((imagecol - start_x2) * (imagecol - start_x2) +
					(imagerow - start_y) * (imagerow - start_y)) <= r * r,
				ifelse(label == 'label2', 'overlap', 'label3'), label))
		labels = vapply(sim$label, function(x) {
			if (x == 'overlap') {
				sample(c('label2', 'label3'), 1)
			} else {
				x
			}
		},
			as.character(1)
		)
		labels = vapply(labels, function(x) {
			if (x %in% c('label2', 'label3')) {
				sample(c('label1', x), 1)
			} else {
				x
			}
		},
			as.character(1)
		)
		sim$label = labels
	} else if (perturbation == 'flag') {
		start_x1 = 254
		start_y1 = 180
		start_x2 = 254
		start_y2 = 180
		start_x3 = 100
		start_y3 = 290
		start_x4 = 100
		start_y4 = 100
		start_x5 = 408
		start_y5 = 290
		start_x6 = 408
		start_y6 = 100
		# r = 50
		# r2 = 25
		r = 80
		r2 = 40

		sim %<>%
			mutate(label = ifelse(
				((imagecol - start_x1) * (imagecol - start_x1) +
					(imagerow - start_y1) * (imagerow - start_y1)) <= r * r,
				'label2', 'background')) %>%
			mutate(label = ifelse(
				((imagecol - start_x2) * (imagecol - start_x2) +
					(imagerow - start_y2) * (imagerow - start_y2)) <= r2 * r2, 'label4', label)) %>%
			mutate(label = ifelse(
				((imagecol - start_x3) * (imagecol - start_x3) +
					(imagerow - start_y3) * (imagerow - start_y3)) <= r2 * r2,
				'label3', label)) %>%
			mutate(label = ifelse(
				((imagecol - start_x4) * (imagecol - start_x4) +
					(imagerow - start_y4) * (imagerow - start_y4)) <= r2 * r2,
				'label2', label)) %>%
			mutate(label = ifelse(
				((imagecol - start_x5) * (imagecol - start_x5) +
					(imagerow - start_y5) * (imagerow - start_y5)) <= r2 * r2,
				'label4', label)) %>%
			mutate(label = ifelse(
				((imagecol - start_x6) * (imagecol - start_x6) +
					(imagerow - start_y6) * (imagerow - start_y6)) <= r2 * r2,
				'label5', label))
		labels = vapply(sim$label, function(x) {
			if (x != 'background') {
				sample(c('label1', x), 1)
			} else {
				x
			}
		},
			as.character(1)
		)
		sim$label = labels
	} else if (perturbation == 'stripes') {
		n_slices = 12
		start_x_vec = seq(range(sim$imagecol)[1], range(sim$imagecol)[2], length.out=n_slices)
		sim$label = 'background'
		label_count = 1
		for (i in c(3, 6, 9)) {
			label_count = label_count + 1
			start_x = start_x_vec[i]
			end_x = start_x_vec[i+1]
			sim %<>%
				mutate(
					label = ifelse(sim$imagecol > start_x & sim$imagecol <= end_x, paste0('label', label_count), label)
				)
		}
		labels = vapply(sim$label, function(x) {
			if (x != 'background') {
				sample(c('label1', x), 1)
			} else {
				x
			}
		},
			as.character(1)
		)
		sim$label = labels
	} 
	# shuffle the points to aid in plotting later
	sim %<>% sample_n(nrow(.))
	return (sim)
}