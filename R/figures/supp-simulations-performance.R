setwd('~/git/vespucci')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(cetcolor)
library(lawstat)
library(nparcomp)
source('R/theme.R')

# read region accuracies ####
reg = readRDS('data/simulations/spatial-acc-summary.rds') %>%
	type_convert() %>%
	filter(
		input != 'circle',
		type == 'AUPRC'
	) %>%
	mutate(
		method = ifelse(method == 'vespucci', 'Meta-learning', 'Exhaustive search')
	)

time = readRDS('data/simulations/time-summary.rds') %>%
	mutate(
			method = ifelse(method == 'vespucci', 'Meta-learning', 'Exhaustive search')
	) %>%
	filter(
			input != 'circle'
	)

other_de_stats = readRDS('data/simulations/summaries/de-results.rds') %>%
    filter(
        input != 'circle'
    )

color_set = data.frame(
	de_method = c(
		'DE only',
		'sparkx',
		'spacgn',
		'cside',
		'moransi',
		'wilcox',
		'nnsvg',
		'spatialDE',
		'heartsvg',
		'squidpy_permutation',
		'squidpy_normality',
		'squidpy_normal_approx_permutation',
		'scran',
		'mast',
		'binSpect_kmeans',
		'binSpect_rank',
		'vespucci'
	),
	color = c(
		'DE only',
		'SPARK-X',
		'SpaGCN',
		'C-SIDE',
		rep('Seurat',2),
		'nnSVG',
		'SpatialDE',
		'HEARTSVG',
		rep('SquidPy',3),
		rep('Giotto', 4),
		'Vespucci'
	),
	x_name = c(
		'NBGMM',
		'SPARK-X',
		'SpaGCN',
		'C-SIDE',
		'Moransi test',
		'Wilcoxon rank-sum test',
		'nnSVG',
		'SpatialDE',
		'HEARTSVG',
		'Squidpy (permutation test)',
		'Squidpy (normality ass.)',
		'Squidpy (normal approx.)',
		'scran',
		'MAST',
		'binSpect (k-means)',
		'binSpect (rank)',
		'Vespucci'
	)
)

# iterate through simulations ####
simulations = c('stripes', 'circle_overlap', 'flag')
for (simulation in simulations) {
	# read simulation object
	if (simulation == 'circle_overlap') {
		sim_file = paste0('data/simulations/objects/input=', simulation,
											'-de_prob=0.5-de_size=2-reps=3-rep_de_jitter=1-rep_depth_jitter=0.3-seed=0.rds')
	} else if (simulation == 'stripes') {
		sim_file = paste0('data/simulations/objects/input=', simulation,
											'-de_prob=0.5-de_size=2-reps=3-rep_de_jitter=1-rep_depth_jitter=0.3-seed=0.rds')
	} else if (simulation == 'flag') {
		sim_file = paste0('data/simulations/objects/input=', simulation,
											'-de_prob=0.5-de_size=2-reps=3-rep_de_jitter=1-rep_depth_jitter=0.3-seed=0.rds')
	}
	sc = readRDS(sim_file)
	meta = sc@meta.data
	
	# read vespucci results
	ves_file = paste0('data/simulations/vespucci/',
										basename(sim_file) %>% gsub("\\.rds", "", .),
										'-ves_seed=42.rds')
	ves_res = readRDS(ves_file)
	spatial_auc_res =  ves_res$spatial_auc_result

	# merge with metadata from sc
	dat0 = spatial_auc_res$aucs %>%
		left_join(meta %>% dplyr::select(barcode, x, y, label_sim))
	dat0 %<>%
		mutate(label_sim = ifelse(label_sim == 'label1', 'background', label_sim))
	
	# set factor levels
	if (simulation == 'circle_overlap') {
		dat0 %<>%
			mutate(label_sim = factor(label_sim, 
				levels = c('background', 'label3', 'label2')))
	} else if (simulation == 'stripes') {
		dat0 %<>%
			mutate(
				label_sim = ifelse(label_sim == 'label4', 'label2', label_sim),
				label_sim = factor(label_sim, levels=c('background', 'label2', 'label3'))
			)
	} else if (simulation == 'flag') {
		dat0 %<>%
			mutate(
				label_sim = ifelse(label_sim == 'label1', 'background', label_sim),
				label_sim = factor(label_sim, levels=c('background', 'label2', 'label3', 'label4', 'label5'))
			)
	}
	
	#############################################################################-
	## a. ground truth ####
	#############################################################################-
	
	truth_pal = pals::kovesi.linear_gow_60_85_c27(100) %>% rev()
	truth_pal = pals::kovesi.linear_blue_95_50_c20(100)
	brks = c(1, 3)
	labels = c('min', 'max')
	p1 = dat0 %>%
		mutate(label_sim = as.integer(label_sim)) %>%
		# arrange(label_sim) %>%
		ggplot(aes(x = y, y = x, fill = label_sim, color = label_sim)) +
		ggtitle('Ground truth') +
		ggrastr::rasterise(
			geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
		) +
		scale_y_continuous(expand = c(0, 0)) +
		scale_x_continuous(expand = c(0, 0)) +
		scale_color_gradientn(colours = truth_pal, name = 'Perturbation   ',
													breaks = brks, labels = labels) +
		scale_fill_gradientn(colours = truth_pal, name = 'Perturbation   ', 
												 breaks = brks, labels = labels) +
		guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE,
																 title.position = "left"),
					 color = guide_colorbar(frame.colour = 'black', ticks = FALSE, 
																	title.position = "left")) +
		coord_fixed() +
		boxed_theme(size_sm = 5, size_lg = 6) +
		theme(axis.title.x = element_blank(),
					axis.title.y = element_blank(),
					axis.text.y = element_blank(),
					axis.text.x = element_blank(),
					axis.ticks.x = element_blank(),
					axis.ticks.y = element_blank(),
					axis.ticks.length.x = unit(0, 'lines'),
					axis.ticks.length.y = unit(0, 'lines'),
					legend.position = 'top',
					legend.justification = 'bottom',
					legend.key.width = unit(0.18, 'lines'),
					legend.key.height = unit(0.18, 'lines')
		)
	p1
	
	#############################################################################-
	## b. Vespucci AUC ####
	#############################################################################-
	
	# interpolate in 2D
	fit = loess(auc ~ x * y, data = dat0, span = 0.015)
	dat0$auc_fit = predict(fit, dat0)
	
	range = range(dat0$auc_fit)
	brks = c(range[1] + 0.1 * diff(range),
					 range[2] - 0.1 * diff(range))
	labels = format(range, digits = 2)
	labels = c(paste0(labels[1], ' '),
						 paste0(' ', labels[2]))
	auc_pal = nr_heat_red_spatial
	auc_pal = pals::kovesi.linear_grey_10_95_c0(100) %>% rev
	auc_pal = pals::kovesi.linear_kryw_5_100_c67(100) %>% rev
	auc_pal = cet_pal(100, name = 'l19') %>% rev()
	
	p2 = dat0 %>%
		# arrange(auc_fit) %>%
		ggplot(aes(x = y, y = x, fill = auc_fit, color = auc_fit)) +
		ggtitle('Vespucci') +
		ggrastr::rasterise(
			geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
		) +
		scale_y_continuous(expand = c(0, 0)) +
		scale_x_continuous(expand = c(0, 0)) +
		scale_fill_gradientn(colours = auc_pal,
												 name = 'AUC   ', labels = labels,
												 limits = range, breaks = brks) +
		scale_color_gradientn(colours = auc_pal,
													name = 'AUC   ', labels = labels,
													limits = range, breaks = brks) +
		guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE),
					 color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
		coord_fixed() +
		boxed_theme(size_lg = 6, size_sm = 5) +
		theme(axis.title.x = element_blank(),
					axis.title.y = element_blank(),
					# axis.text.y = element_text(angle = 90, hjust = c(1, 0)),
					# axis.text.x = element_text(hjust = c(0, 1)),
					axis.text.y = element_blank(),
					axis.text.x = element_blank(),
					axis.ticks.x = element_blank(),
					axis.ticks.y = element_blank(),
					axis.ticks.length.x = unit(0, 'lines'),
					axis.ticks.length.y = unit(0, 'lines'),
					legend.position = 'top',
					# legend.justification = 'bottom',
					legend.key.width = unit(0.18, 'lines'),
					legend.key.height = unit(0.18, 'lines'),
					plot.title = element_text(size = 6))
	p2
	
	#############################################################################-
	## c. Vespucci vs. Magellan, walltime ####
	#############################################################################-
	
	time0 = time %>% 
		filter(input == simulation)
	labs = time0 %>%
		group_by(input, method) %>%
		summarise(
			time = max(time),
			label = paste0(round(mean(time), 1), ' h')
		)
	
	# paired t-test
	delta = time0 %>% 
		group_by(seed) %>% 
		summarise(delta = time[method == 'Meta-learning'] - 
								time[method == 'Exhaustive search']) %>% 
		ungroup()
	pval = t.test(delta$delta)$p.value %>% 
		format(format = 'f', digits = 2) %>% 
		paste0('p = ', .)
	pval_df = data.frame(label = pval)
	
	range = boxplot(time ~ method, dat = time0)$stats %>% range
	og_pal = c('Exhaustive search' = nr_base_5[2],
						 'Meta-learning' = nr_base_5[1])
	p3 = time0 %>%
		ggplot(aes(x = method, y = time)) +
		# facet_grid(~ input) +
		geom_label(data = pval_df, aes(x = 1.5, y = -Inf, label = '____________________________'),
							 size = 1.75, vjust = 0, label.size = NA, color = NA, fill = 'grey96',
							 label.padding = unit(0.2, 'lines')) +
		geom_label(data = pval_df, aes(x = 1.5, y = -Inf, label = label),
							 size = 1.5, vjust = 0, label.size = NA, fill = NA,
							 label.padding = unit(0.45, 'lines')) +
		geom_boxplot(aes(color = method, fill = method), 
								 alpha = 0.4, width = 0.6, size = 0.35, outlier.shape = NA) +
		# geom_jitter(shape = 21, size = 0.4, stroke = 0.25, height = 0, 
		#             width = 0.15) +
		geom_text(dat = labs, aes(y = time, x = method, label = label),
							color = 'black',
							size = 1.5,
							vjust=-1) +
		boxed_theme() +
		scale_color_manual(values = og_pal) +
		scale_fill_manual(values = og_pal) +
		scale_x_discrete(labels = ~ gsub(" ", '\n', .)) +
		scale_y_continuous('Runtime, hours', breaks = pretty_breaks(4), 
											 expand = expansion(c(0.25, 0.25))) +
		coord_cartesian(ylim = range) +
		theme(
			aspect.ratio = 1.7,
			axis.title.x = element_blank(),
			axis.text.x = element_text(angle = 45, hjust = 1, lineheight = 0.8),
			legend.position = 'none',
			legend.justification = 'bottom',
			legend.title = element_blank()) 
	p3
	
	#############################################################################-
	## d. Vespucci vs. Magellan, region AUC ####
	#############################################################################-
	
	reg0 = reg %>% 
		filter(input == simulation)
	labs = reg0 %>%
		group_by(method, input) %>%
		summarise(
			val = max(val),
			label = round(mean(val), 3)
		)
	
	# paired t-test
	delta = reg0 %>% 
		group_by(seed) %>% 
		summarise(delta = val[method == 'Meta-learning'] - 
								val[method == 'Exhaustive search']) %>% 
		ungroup()
	pval = t.test(delta$delta)$p.value %>% 
		formatC(format = 'f', digits = 2) %>% 
		paste0('p = ', .)
	pval_df = data.frame(label = pval)
	
	range = boxplot(val ~ method, dat = reg0)$stats %>% range
	og_pal = c('Exhaustive search' = nr_base_5[2],
						 'Meta-learning' = nr_base_5[1])
	p4 = reg0 %>%
		ggplot(aes(x = method, y = val)) +
		# facet_grid(~ input) +
		geom_label(data = pval_df, aes(x = 1.5, y = -Inf, label = '____________________________'),
							 size = 1.75, vjust = 0, label.size = NA, color = NA, fill = 'grey96',
							 label.padding = unit(0.2, 'lines')) +
		geom_label(data = pval_df, aes(x = 1.5, y = -Inf, label = label),
							 size = 1.5, vjust = 0, label.size = NA, fill = NA,
							 label.padding = unit(0.45, 'lines')) +
		geom_boxplot(aes(color = method, fill = method), 
								 alpha = 0.4, width = 0.6, size = 0.35, outlier.shape = NA) +
		# geom_jitter(shape = 21, size = 0.4, stroke = 0.25, height = 0, 
		#             width = 0.15) +
		geom_text(dat = labs, aes(y = val, x = method, label = label),
							color = 'black',
							size = 1.5,
							vjust=-1) +
		boxed_theme() +
		scale_color_manual(values = og_pal) +
		scale_fill_manual(values = og_pal) +
		scale_x_discrete(labels = ~ gsub(" ", '\n', .)) +
		scale_y_continuous('AUPRC', breaks = pretty_breaks(4), 
											 expand = expansion(c(0.25, 0.25))) +
		coord_cartesian(ylim = range) +
		theme(
			aspect.ratio = 1.7,
			axis.title.x = element_blank(),
			axis.text.x = element_text(angle = 45, hjust = 1, lineheight = 0.8),
			legend.position = 'none',
			legend.justification = 'bottom',
			legend.title = element_blank()) 
	p4
	
	#############################################################################-
	## e. AUPR boxplot ####
	#############################################################################-
	ves_de_res = ves_res$de_feature_result

	DE = rbind(
		ves_de_res,
		other_de_stats %>% filter(input == !!input)
	) %>% dplyr::rename(val = auprc_integral)
	DE %>% dplyr::select(seed, de_method) %>% table()
	DE %<>% left_join(color_set, by = 'de_method')

	# test differences ####
	pairs = tidyr::crossing(method1 = unique(DE$de_method), method2 = unique(DE$de_method)) %>% 
		filter(method1 != method2)

	tests = pmap_dfr(pairs, function(...) {
		current = tibble(...)
		print(current)
		
		vec1 = filter(DE, de_method == current$method1) %>% 
			arrange(input, seed) %>% pull(val)
		vec2 = filter(DE, de_method == current$method2) %>% 
			arrange(input, seed) %>% pull(val)
		t = t.test(vec1, vec2)$p.value
		pt = t.test(vec1, vec2, paired = TRUE)$p.value
		w = wilcox.test(vec1, vec2)$p.value
		pw = wilcox.test(vec1, vec2, paired = TRUE)$p.value
		bm = brunner.munzel.test(vec1, vec2)$p.value
		pbm_df = data.frame(method = rep(c('method1', 'method2'), each = length(vec1)),
												idx = rep(seq_along(vec1), 2),
												value = c(vec1, vec2)) %>% 
			arrange(idx)
		pbm = npar.t.test.paired(value ~ method, pbm_df, nperm = 0)$Analysis['BM', 'p.value']
		data.frame(test = c('t-test', 'paired t-test', 'wilcox', 'paired wilcox',
												'brunner-munzel', 'paired brunner-munzel'),
							 pval = c(t, pt, w, pw, bm, pbm)) %>% 
			cbind(current, .)
	})

	# calculate delta-mean and median
	deltas = pmap_dfr(pairs, function(...) {
		current = tibble(...)
		vec1 = filter(DE, de_method == current$method1) %>% 
			arrange(seed) %>% pull(val)
		vec2 = filter(DE, de_method == current$method2) %>% 
			arrange(seed) %>% pull(val)
		median = median(vec2 - vec1)
		mean = mean(vec2 - vec1)
		mutate(current, delta_median = median, delta_mean = mean)
	})
	
	# set color palette
	pal = nr_base_11_light[1:10]
	names(pal) = unique(DE$color[!DE0$color %in% c('Vespucci')])
	pal['Vespucci'] = nr_base_5[1]
	alpha_pal = c('max'=1, 'min'=0.3)
	
	# add OOT methods
	oot_df = data.frame(x_name = c('SPARK', 'SPADE', 'trendsceek', 'GPCounts'), input = DE$input %>% first) %>% 
		mutate(color = x_name) %>% 
		filter(!x_name %in% DE$x_name)
	DE %<>% bind_rows(oot_df)
	
	labs = DE %>%
		group_by(x_name, input, color) %>%
		summarize(
			# stats_val = median(val),
			stats_val = mean(val),
			val = max(val)
		) %>%
		ungroup() %>%
		mutate(
			text_val = ifelse(stats_val < 0, 'NA', format(stats_val, digits = 2)),
			text_y = ifelse(val < 0, 0, val)
		) %>% 
		replace_na(list(text_val = 'OOT', text_y = -Inf))
		
	med = function(x) stats::median(x) %>% replace(is.na(.), 0)
	p5 = DE %>%
		ggplot(aes(x = reorder(x_name, val, med), y = val,
							 fill = color, color = color)) +
		geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA) +
		geom_label(dat = labs,
							 aes(label = text_val, y = text_y), color = 'black',
							 label.padding = unit(0.35, 'lines'),
							 label.size = NA, fill = NA,
							 size = 1.75, hjust = 0, vjust = 0.5,
							 show.legend = FALSE) +
		coord_flip() +
		scale_y_continuous('AUPR', breaks = pretty_breaks(),
											 expand = expansion(c(0.03, 0.125))) +
		scale_color_manual('', values = pal) +
		scale_fill_manual('', values = pal) +
		boxed_theme() +
		theme(axis.title.y = element_blank(),
					aspect.ratio = 1.5,
					legend.position = 'none',
					legend.justification = 'bottom',
					legend.key.width = unit(0.18, 'lines'),
					legend.key.height = unit(0.15, 'lines')
		)
	p5
	
	#############################################################################-
	## f. delta-AUPR heatmap ####
	#############################################################################-
	
	lvls = DE %>% 
		drop_na(val) %$%
		reorder(de_method, val, mean) %>% 
		levels()
	range = range(deltas$delta_mean)
	brks = c(range[1] + 0.1 * diff(range),
					 range[2] - 0.1 * diff(range))
	xlab = with(color_set, setNames(x_name, de_method))
	labels = tests %>% 
		filter(test == 'paired t-test') %>% 
		mutate(method1 = factor(method1, levels = lvls),
					 method2 = factor(method2, levels = lvls),
					 lab = ifelse(pval < 0.001, '***',
												ifelse(pval < 0.01, '**',
															 ifelse(pval < 0.05, '*', ''))))
	p6 = deltas %>% 
		mutate(method1 = factor(method1, levels = lvls),
					 method2 = factor(method2, levels = lvls)) %>% 
		ggplot(aes(x = method1, y = method2)) +
		# geom_tile(color = 'white', aes(fill = delta_median)) +
		geom_tile(color = 'white', aes(fill = delta_mean)) +
		geom_text(data = labels, size = 1.5, aes(label = lab), nudge_y = -0.15) +
		scale_x_discrete(expand = c(0, 0), labels = xlab) +
		scale_y_discrete(expand = c(0, 0), labels = xlab) +
		scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
													 name = expression(Delta~AUPR),
													 breaks = brks,
													 labels = format(range, digits = 2)) +
		guides(fill = guide_colorbar(ticks = FALSE, frame.colour = 'black')) +
		coord_fixed() +
		boxed_theme() +
		theme(axis.title.x = element_blank(),
					axis.title.y = element_blank(),
					# axis.text.x = element_blank(),
					# axis.ticks.x = element_blank(),
					axis.text.x = element_text(angle = 45, hjust = 1),
					legend.key.width = unit(0.18, 'lines'),
					legend.key.height = unit(0.15, 'lines'),
					legend.position = 'right',
					legend.justification = 'bottom')
	p6
	
	#############################################################################-
	## combine and save ####
	#############################################################################-
	
	# combine
	row1 = p1 + p2 + p3 + p4 + plot_layout(nrow = 1)
	ggsave(paste0("fig/final/EFig2/", simulation, "-row1.pdf"), row1,
				 width = 12, height = 5, units = "cm", useDingbats = FALSE)
	ggsave(paste0("fig/final/EFig2/", simulation, "-row2.pdf"), p5,
				 width = 8, height = 6.3, units = "cm", useDingbats = FALSE)
	ggsave(paste0("fig/final/EFig2/", simulation, "-row3.pdf"), p6,
				 width = 9, height = 6.6, units = "cm", useDingbats = FALSE)
}

