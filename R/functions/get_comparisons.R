## Get subset expression matrices containing all comparisons for a given
## dataset.
get_comparisons = function(dataset, expr, meta) {
	# set up container
	results = list()
	if (dataset == 'Squair2023_3d') {
		message('  processing comparison: ', 'uninjured vs 7days', '...')
		meta0 = meta %>%
			mutate(idx = row_number(), label = group) %>%
			filter(group %in% c("uninjured", "7days"))
		expr0 = expr %>% extract(, meta0$idx)
		results[['uninjured-7days']] = list(expr = expr0, meta = meta0)
	} 
	else if (dataset == 'Maniatis2019') {
		timepoints = c('p100')
		for (timepoint in timepoints) {
			key = paste0("WT-SOD|", timepoint)
			message('  processing comparison: ', key, ' ...')
			meta0 = meta %>% 
				mutate(idx = row_number()) %>%
				filter(label %in% c("B6SJLSOD1-G93A", "B6SJLSOD1-WT")) %>%
				filter(timepoint == !!timepoint)
			expr0 = expr %>% extract(, meta0$idx)
			results[[key]] = list(expr = expr0, meta = meta0)  
		}
	}
	else if (dataset == 'Kathe2022') {
		# only one comparison
		grid = data.frame('group1'='EES_REHAB', 'group2'='SCI')

		for (grid_idx in seq_len(nrow(grid))) {
			group1 = grid$group1[grid_idx]
			group2 = grid$group2[grid_idx]
			key = paste0(group1, '-', group2)
			message('  processing comparison: ', key, ' ...')

			meta0 = meta %>%
				mutate(idx = row_number()) %>%
				filter(label %in% c(group1, group2))
				expr0 = expr %>% extract(, meta0$idx)
				results[[key]] = list(expr = expr0, meta = meta0)
		}
	} 
	else if (dataset == 'Zeng2023') {
		message('  processing comparison: ', 'disease vs control at 13 months', '...')
		meta0 = meta %>%
			mutate(
				idx = row_number(),
				label = group
			) %>%
			filter(time == '13months')
		expr0 = expr %>% extract(, meta0$idx)
		results[['Disease-Control|13months']] = list(expr = expr0, meta = meta0)
	} 
	else if (dataset == 'Calcagno2022') {
		meta %<>%
			mutate(
				label = timepoint
			)
		labels = c('d1', 'd7')
		comparison_grid = tidyr::crossing(label1 = labels, label2 = labels)

		comparison_grid %<>% 
			filter(label2 > label1)

		for (i in 1:nrow(comparison_grid)) {
			tmp_row = comparison_grid[i,]
			comparison_name = paste0(tmp_row$label1, '-', tmp_row$label2)
			message('  processing comparison: ', paste0(tmp_row$label1, ' vs ', tmp_row$label2), '...')
			meta0 = meta %>%
				mutate(
					idx = row_number()
				) %>%
				filter(label %in% c(tmp_row$label1, tmp_row$label2))
			expr0 = expr %>% extract(, meta0$idx)
			results[[comparison_name]] = list(expr = expr0, meta = meta0)
		}
	} else {
		stop("invalid dataset: ", dataset, " ...")
		}
		# drop all unused factor levels
		for (comparison_idx in seq_along(results)) {
			results[[comparison_idx]]$meta %<>% droplevels()
	}
	return(results)
}