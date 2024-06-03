
### all modified functions

splat_simulate = function(
	params = newSplatParams(),
	preset_sim_metadata = NULL,
	multi_de_list = NULL,
	de_specific_cells_list = NULL,
	group.prob = c(0.5, 0.5),
	method = "groups",
	batchCells = 10000,
	n_reps_per_group = 3,
	rep_depth_jitter = 1,
	seed = 42,
	false_DE = FALSE,
	verbose = TRUE
	) {
	# Set random seed
	set.seed(seed)

	# Get the parameters we are going to use
	nCells = getParam(params, "nCells")
	nGenes = getParam(params, "nGenes")
	nBatches = getParam(params, "nBatches")
	batch.cells = getParam(params, "batchCells")
	nGroups = getParam(params, "nGroups")
	group.prob = getParam(params, "group.prob")

	if (nGroups == 1 && method == "groups") {
		warning("nGroups is 1, switching to single mode")
		method = "single"
	}

	# Set up name vectors
	if (verbose) {message("Creating simulation object...")}
	cell.names = paste0("Cell", seq_len(nCells))
	gene.names = paste0("Gene", seq_len(nGenes))
	batch.names = paste0("Batch", seq_len(nBatches))
	if (method == "groups") {
		group.names = paste0("Group", seq_len(nGroups))
	} else if (method == "paths") {
		group.names = paste0("Path", seq_len(nGroups))
	}

	# Create SingleCellExperiment to store simulation
	cells =  data.frame(Cell = cell.names)
	rownames(cells) = cell.names
	features = data.frame(Gene = gene.names)
	rownames(features) = gene.names
	sim = SingleCellExperiment(rowData = features, colData = cells, metadata = list(Params = params))

	if (!is.null(preset_sim_metadata)){
		rownames(preset_sim_metadata) = preset_sim_metadata$Cell
		preset_sim_metadata = preset_sim_metadata[colData(sim)$Cell,]
		for (colname in colnames(preset_sim_metadata)){
			colData(sim)@listData[[colname]] = preset_sim_metadata[,colname]
		}
	}

	# Set replicates
	if (n_reps_per_group > 1){
		colData(sim)@listData[["replicate"]] = paste0('Rep', sample(1:n_reps_per_group, ncol(sim), replace=T))
	} else {
		colData(sim)@listData[["replicate"]] = 'Rep1'
	}

	# Make batches vector which is the index of param$batchCells repeated
	# params$batchCells[index] times
	batches = lapply(seq_len(nBatches), function(i, b) {rep(i, b[i])},
					b = batch.cells)
	batches = unlist(batches)
	colData(sim)$Batch = batch.names[batches]

	## only if we decide the metadata already
	if (is.null(multi_de_list) & is.null(preset_sim_metadata)){
		if (method != "single") {
			groups = sample(seq_len(nGroups), nCells, prob = group.prob,
							replace = TRUE)
			colData(sim)$Group = factor(group.names[groups], levels = group.names)
		}	
	}

	if (verbose) {message("Simulating library sizes...")}
		sim = splatter:::splatSimLibSizes(sim, params)
	if (verbose) {message("Simulating gene means...")}
	sim = splatter:::splatSimGeneMeans(sim, params)
	if (nBatches > 1) {
		if (verbose) {message("Simulating batch effects...")}
		sim = splatSimBatchEffects(sim, params)
	}
	sim = splatter:::splatSimBatchCellMeans(sim, params)

	# custom set for now
	all_genes = rownames(sim)
	n_genes = nrow(sim)
	de_genes = sample(all_genes, floor(n_genes*multi_de_list[['Perturbation']]$de_prob))
	
	n_reps_per_group_2 = ifelse(false_DE, 1, n_reps_per_group)

	if (method == "single") {
		sim = splatSimSingleCellMeans(sim, params)
	} else if (method == "groups") {
		## original DE effect
		for (de_label_colname in names(multi_de_list)){
			sim = splatSimGroupDE(sim, params, genes = de_genes, de_label_colname = de_label_colname, custom_de_details = multi_de_list[[de_label_colname]], n_reps_per_group = n_reps_per_group_2, verbose = verbose)
		}
		if (verbose) {message("Simulating cell means...")}
		sim = splatSimGroupCellMeans(sim, params, multi_de_list = multi_de_list, de_specific_cells_list = de_specific_cells_list, n_reps_per_group = n_reps_per_group, verbose = verbose)
	} else {
		if (verbose) {message("Simulating path endpoints...")}
		sim = splatSimPathDE(sim, params)
		if (verbose) {message("Simulating path steps...")}
		sim = splatSimPathCellMeans(sim, params)
	}

	if (verbose) {message("Simulating BCV...")}
	sim = splatter:::splatSimBCVMeans(sim, params)
	if (verbose) {message("Simulating counts...")}
	sim = splatSimTrueCounts(sim, params, rep_depth_jitter=rep_depth_jitter, n_reps_per_group = n_reps_per_group)
	if (verbose) {message("Simulating dropout (if needed)...")}
	sim = splatter:::splatSimDropout(sim, params)

	sc = sim %>% logNormCounts() %>% as.Seurat

	sc@meta.data %<>%
		left_join(sim_grid %>% 
			dplyr::rename(x = imagerow, y = imagecol) %>%
			dplyr::select(x,y,Cell), 
		by='Cell')

	if (!false_DE) {
		set_perturb_labels = vapply(sim_grid$label, function(x) {
			ifelse(x == 'background', 
				sample(c('Perturbation_1', 'Perturbation_2'), 1),
				gsub('label', 'Perturbation_', x)
			)
			},
			as.character(1)
		)
		set_perturb_labels = unname(set_perturb_labels)
		set_perturb_labels[set_perturb_labels != 'Perturbation_1'] = 'Perturbation_2'

		rownames(sc@meta.data) = sc@meta.data$Cell
		sc@meta.data = sc@meta.data[colnames(sc),]
		sc@meta.data %<>%
			mutate(
				barcode = Cell,
				label_sim = label,
				label = set_perturb_labels,
				replicate = paste0(gsub('_', '', as.character(label)), '_', replicate)
		)
	} else {
		rownames(sc@meta.data) = sc@meta.data$Cell
		sc@meta.data = sc@meta.data[colnames(sc),]
		sc@meta.data %<>%
			mutate(
				ori_rep = replicate,
				label_sim = Perturbation,
				replicate = gsub('Perturbation', 'Rep', as.character(Perturbation)),
				barcode = Cell
			)

		control_reps = sample(unique(sc@meta.data$replicate), n_reps_per_group)
		sc@meta.data %<>%
			mutate(label = ifelse(replicate %in% control_reps, 'Perturbation_1', 'Perturbation_2'))
		sc@meta.data$label = factor(sc@meta.data$label)
	}
	return(sc)
}

splatSimGroupDE = function(sim, params, genes = NULL, de_label_colname = NULL, custom_de_details = NULL, n_reps_per_group = 1, verbose = FALSE) {

	nGenes = getParam(params, "nGenes")
	nGroups = getParam(params, "nGroups")
	de.prob = getParam(params, "de.prob")
	de.downProb = getParam(params, "de.downProb")
	de.facLoc = getParam(params, "de.facLoc")
	de.facScale = getParam(params, "de.facScale")
	means.gene = rowData(sim)$GeneMean

	if (is.null(custom_de_details)){
		for (idx in seq_len(nGroups)) {
			de.facs = getLNormFactors(nGenes, de.prob[idx], de.downProb[idx], de.facLoc[idx], de.facScale[idx])
			group.means.gene = means.gene * de.facs
			rowData(sim)[[paste0("DEFacGroup", idx)]] = de.facs
		}	
	} else {
		de_prob = custom_de_details[['de_prob']]
		de_effect_size_vec = custom_de_details[['de_effect_size']]
		de_fac_scale = custom_de_details[['de_fac_scale']]
		de_down_prob = custom_de_details[['de_down_prob']]
		selected_de_idxs = custom_de_details[['de_selected_idx']]
		de_direction = custom_de_details[['de_direction']]
		rep_de_facloc_jitter_vec = custom_de_details[['rep_de_facloc_jitter']]

		if (is.factor(colData(sim)[[de_label_colname]])) {
			de_labels = levels(colData(sim)[[de_label_colname]])
		} else {
			de_labels = sort(unique(colData(sim)[[de_label_colname]]))	
		}

		if (length(de_effect_size_vec) == 1) {
			de_effect_size_vec = rep(de_effect_size_vec, length(de_labels))
		}

		if (length(rep_de_facloc_jitter_vec) == 1) {
			rep_de_facloc_jitter_vec = rep(rep_de_facloc_jitter_vec, length(de_labels))
		}

		names(de_effect_size_vec) = de_labels
		names(rep_de_facloc_jitter_vec) = de_labels

		selected_de_labels = as.character(de_labels[selected_de_idxs])
		
		nGroups_for_genes = ifelse(is.null(selected_de_idxs), length(de_labels), length(selected_de_idxs))
		
		if (is.null(genes)) genes = sample(rownames(sim), floor(de_prob * nrow(sim)))

		if (nGroups_for_genes > 1 & length(genes) > 0){
			genes_by_group = map(1:nGroups_for_genes, function(x){genes})
		} else {
			genes_by_group = list(genes)
		}
		if (!is.null(selected_de_idxs) & length(genes) > 0){
			names(genes_by_group) = selected_de_labels
		}

		if (!is.null(rowData(sim)$DEGeneMean)){
			means.gene = rowData(sim)$DEGeneMean
		}

		same_genes_direction = NULL
		for (de_label in de_labels){
			if (de_label %in% selected_de_labels){
				
				de_effect_size  = de_effect_size_vec[de_label]
				rep_de_facloc_jitter = rep_de_facloc_jitter_vec[de_label]

				genes_to_run = rownames(sim) %in% genes_by_group[[de_label]] 

				genes_direction = NULL
				if (!is.null(de_direction)){
					if (de_direction == 'up') {
						genes_direction = rep(1, sum(genes_to_run))
					} else if (de_direction == 'down'){
						genes_direction = rep(-1, sum(genes_to_run))
					}
				}

				if (n_reps_per_group == 1){

					if (verbose) message(paste0('Adding de effect to ', de_label, ' with effect size ', de_effect_size, ' at prob ', de_prob))
					
					de.facs = getLNormFactors(nGenes, de_prob, de_down_prob, de_effect_size, de_fac_scale, genes = genes_to_run, genes_direction = genes_direction)
					# group.means.gene = means.gene * de.facs

					rowData(sim)[[de_label]] = de.facs$factors
					rowData(sim)[[paste0(de_label, "_is_selected")]] = de.facs$is_selected
					rowData(sim)[[paste0(de_label, "_direction")]] = de.facs$direction
				} else {
					if (is.null(rep_de_facloc_jitter)){
						rep_de_facloc_jitter = 0
					}
					n_reps_per_group = length(unique(colData(sim)$replicate))

					genes_direction = NULL
					de_facLoc_vector = rep(de_effect_size, n_reps_per_group)
					de_facLoc_jitter = seq_len(n_reps_per_group) * rep_de_facloc_jitter
					de_facLoc_jitter = de_facLoc_jitter - median(de_facLoc_jitter)
					de_facLoc_vector = de_facLoc_vector + de_facLoc_jitter
					# minimum is zero
					de_facLoc_vector[de_facLoc_vector < 0] = 0
					genes_direction = same_genes_direction

					for (rep_idx in seq_len(n_reps_per_group)) {
						
						rep_de_facs = getLNormFactors(nGenes, de_prob, de_down_prob, de_facLoc_vector[rep_idx], de_fac_scale, genes = genes_to_run, genes_direction = genes_direction)

						if (rep_idx == 1) {
							genes_direction = rep_de_facs$direction[genes_to_run]
							genes_direction = ifelse(genes_direction == 'up', 1, -1)
							if (is.null(same_genes_direction)) {
									same_genes_direction = genes_direction
								}
						}

						curr_de_label = paste0(de_label, '_Rep', rep_idx)

						if (verbose) message(paste0('Adding de effect to ', curr_de_label, ' with effect size ', de_facLoc_vector[rep_idx], ' at prob ', de_prob))

						rowData(sim)[[curr_de_label]] = rep_de_facs$factors
						rowData(sim)[[paste0(curr_de_label, "_is_selected")]] = rep_de_facs$is_selected
						rowData(sim)[[paste0(curr_de_label, "_direction")]] = rep_de_facs$direction
					}	
				}
				
			} else {
				if (n_reps_per_group <= 1){
					if (verbose) message(paste0('No de effect to ', de_label))
					rowData(sim)[[de_label]] = 1
					rowData(sim)[[paste0(de_label, "_is_selected")]] = FALSE
					rowData(sim)[[paste0(de_label, "_direction")]] = 'n.s.'	
				} else {
					for (rep_idx in 1:n_reps_per_group){
						curr_de_label = paste0(de_label, '_Rep', rep_idx)
						if (verbose) message(paste0('No de effect to ',curr_de_label))
						rowData(sim)[[curr_de_label]] = 1
						rowData(sim)[[paste0(curr_de_label, "_is_selected")]] = FALSE
						rowData(sim)[[paste0(curr_de_label, "_direction")]] = 'n.s.'	
					}
				}
				
			}
		}

	}

	return(sim)
}

getLNormFactors = function(n.facs, sel.prob, neg.prob, fac.loc, fac.scale, genes = NULL, genes_direction = NULL) {

	if (!is.null(genes)) {
		is.selected = genes
	} else {
		is.selected = as.logical(rbinom(n.facs, 1, sel.prob))
	}
	n.selected = sum(is.selected)
	if (is.null(genes_direction)) {
		dir.selected = (-1) ^ rbinom(n.selected, 1, neg.prob)
	} else {
		dir.selected = genes_direction
	}

	facs.selected = rlnorm(n.selected, fac.loc, fac.scale)
	# Reverse directions for factors that are less than one
	dir.selected[facs.selected < 1] = -1 * dir.selected[facs.selected < 1]
	factors = rep(1, n.facs)
	factors[is.selected] = facs.selected ^ dir.selected

	direction = as.character(
		ifelse(factors > 1, 'up', ifelse(factors < 1, 'down', 'n.s.')))

	## we will return a bit more information for the sake of benchmarking
	## NOTE: I will fix the DE / batch functions but changing this will
	## likely break some other functionality we aren't using.
	out = data.frame(
		is_selected = is.selected,
		direction = direction,
		factors = factors
	)

	return(out)
}

splatSimGroupCellMeans = function(sim, params, multi_de_list = NULL, de_specific_cells_list = NULL, n_reps_per_group = 1, verbose = FALSE) {

	nGroups = getParam(params, "nGroups")
	cell.names = colData(sim)$Cell
	gene.names = rowData(sim)$Gene
	groups = colData(sim)$Group
	group.names = levels(groups)
	exp.lib.sizes = colData(sim)$ExpLibSize
	batch.means.cell = assays(sim)$BatchCellMeans

	if (is.null(multi_de_list)){
		group.facs.gene = rowData(sim)[, paste0("DEFac", group.names)]
		cell.facs.gene = as.matrix(group.facs.gene[, paste0("DEFac", groups)])
		cell.means.gene = batch.means.cell * cell.facs.gene
		cell.props.gene = t(t(cell.means.gene) / colSums(cell.means.gene))
		base.means.cell = t(t(cell.props.gene) * exp.lib.sizes)			
	} else {
		out.means.gene = batch.means.cell
		# here we go with the mind-melting matrices
		for (de_label_name in names(multi_de_list)){
			if (verbose) message(paste0('Adding DE facs of ', de_label_name))

			if (n_reps_per_group <= 1){
				cell_level_de_labels = as.character(colData(sim)[[de_label_name]])
				de_label_cols = unique(cell_level_de_labels)
			} else {
				cell_level_de_labels = paste0(as.character(colData(sim)[[de_label_name]]), '_', as.character(colData(sim)[['replicate']]))
				de_label_cols = unique(cell_level_de_labels)
			}

			de_label_cols = de_label_cols[de_label_cols %in% colnames(rowData(sim))]

			de.facs.genes = rowData(sim)[, de_label_cols]
			cell.facs.gene = as.matrix(de.facs.genes[,cell_level_de_labels])

			# remove effect off all non-chosen cells
			if (!is.null(de_specific_cells_list)){
				de_specific_cells = de_specific_cells_list[[de_label_name]]
				de_specific_cells_bool = !cell.names %in% de_specific_cells
				if (verbose) message(paste0('Filtering DE facs of ', sum(de_specific_cells_bool), ' cells out of ', de_label_name))
				cell.facs.gene[, de_specific_cells_bool] = 1
			}

			out.means.gene = out.means.gene * cell.facs.gene
		}
		# cell.props.gene = t(t(out.means.gene) / colSums(out.means.gene))
		# divides by the number of counts 
		cell.props.gene = t(t(out.means.gene) / colSums(batch.means.cell))
		base.means.cell = t(t(cell.props.gene) * exp.lib.sizes)
	}


	colnames(base.means.cell) = cell.names
	rownames(base.means.cell) = gene.names
	assays(sim)$BaseCellMeans = base.means.cell

	return(sim)
}

splatSimTrueCounts = function(sim, params,
	n_reps_per_group = 1,
	rep_depth_jitter = 0
	) {

	cell.names = colData(sim)$Cell
	gene.names = rowData(sim)$Gene
	nGenes = getParam(params, "nGenes")
	nCells = getParam(params, "nCells")
	cell.means = assays(sim)$CellMeans

	true.counts = matrix(rpois(
		as.numeric(nGenes) * as.numeric(nCells),
		lambda = cell.means),
	nrow = nGenes, ncol = nCells)

	colnames(true.counts) = cell.names
	rownames(true.counts) = gene.names

	# add replicate depth jitter here
	if (rep_depth_jitter > 0 & n_reps_per_group > 1) {
		reps = paste0('Rep', 1:n_reps_per_group)
		cell_reps = colData(sim)$replicate

		## downsampleMatrix only go downwards (obviously)
		
		rep.depth.vector = 1 - (rep(rep_depth_jitter, n_reps_per_group) * seq_len(n_reps_per_group)) + rep_depth_jitter
		names(rep.depth.vector) = unique(reps)

		# minimum is 0.01 proportion 
		reps[reps < 0.01] = 0.01
		
		for (rep in reps){
			message(paste0('Downsampling sequencing depth of ', rep, ' to ', rep.depth.vector[rep]*100, '%'))
			reps_bool = cell_reps == rep
			true.counts[, reps_bool] = as.matrix(scuttle::downsampleMatrix(true.counts[, reps_bool], rep.depth.vector[rep]))
		}
	}

	assays(sim)$TrueCounts = true.counts

	return(sim)
}