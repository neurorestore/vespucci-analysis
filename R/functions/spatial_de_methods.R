library(Seurat)
library(tidyverse)
library(magrittr)
library(spacexr)
library(nnSVG)
library(trendsceek) # remotes::install_version('spatstat', version='1.64.1')
library(SpatialExperiment)
library(SPARK)
library(reticulate)
library(SPADE)
use_condaenv('de-methods')
source_python('python/analysis/spade_naiveDE_functions.py')
library(Giotto)

# methods that has filtering: seurat methods, nnsvg

run_spatial_de = function(input, meta, gene_features, de_method = 'wilcox', filter_genes = FALSE, selected_features = NULL, giotto_python_path = NULL){
	if (de_method %in% c('wilcox', 'markvariogram', 'moransi')){
		sc = CreateSeuratObject(input, assay='Spatial')
		sc@meta.data = meta
		coords_df = meta %>%
			select(x,y) %>%
			set_rownames(meta$barcode)
		sc@images$image = new(Class='SlideSeq', assay='Spatial', key='image_', coordinates = coords_df)

		if (de_method == 'wilcox'){
			sc %<>% NormalizeData()
			Idents(sc) = sc@meta.data$label
			labels = unique(meta$label)
			output_res = Seurat::FindMarkers(sc, 
				ident.1 = labels[1], ident.2 = labels[2], 
				test.use='wilcox',
				logfc.threshold = 0,
				min.pct = 0
				)
			output_res %<>%
				rownames_to_column('gene') %>%
				dplyr::select(gene, avg_log2FC, p_val) %>%
				as_tibble()
			output_res$p_val_adj = p.adjust(output_res$p_val, method = 'BH')
		} else {
			features = rownames(sc)
			sc %<>% SCTransform(., assay = "Spatial", variable.features.n = length(features), verbose = TRUE)
			sc %<>% FindSpatiallyVariableFeatures(., assay = "SCT", features = features, selection.method = 'moransi')
						
			if (de_method == 'moransi') {
				output_res = sc$SCT@meta.features
				output_res %<>%
					as_tibble() %>%
					dplyr::select(2, 1) %>%
					set_colnames(c('stat', 'p_val'))
				output_res$gene = rownames(sc$SCT@meta.features)
				output_res %<>%
					relocate(gene, .before = stat) %>%
					filter(!is.na(p_val))
				output_res$p_val_adj = p.adjust(output_res$p_val, method = 'BH')
			} else {
				output_res = sc$SCT@meta.features
				output_res$gene = rownames(sc$SCT@meta.features)
				output_res %<>%
					as_tibble() %>%
					dplyr::select(4, 1) %>%
					set_colnames(c('gene', 'stat')) %>%
					filter(!is.na(stat))
			}
		}
	} else if (de_method == 'nnsvg') {
		## https://github.com/lmweber/nnSVG
		spe = SpatialExperiment(
			assay = input,
			colData = meta,
			spatialCoordsNames = c("x", "y")
		)
		assayNames(spe) = 'counts'
		rowData(spe)$gene_name = rownames(spe)
		if (filter_genes) {
			spe = filter_genes(spe)
			spe = scuttle::computeLibraryFactors(spe)
		}
		spe = scuttle::logNormCounts(spe)
		spe = nnSVG(spe, n_threads = 1)
		output_res = rowData(spe)@listData %>%
			data.frame() %>%
			select(LR_stat, pval, padj) %>%
			set_colnames(c('stat', 'p_val', 'p_val_adj'))
		output_res$gene = rownames(spe)
		output_res %<>%
			relocate(gene, .before=stat)
	} else if (de_method == 'cside') {
		## https://raw.githack.com/dmcable/spacexr/master/vignettes/differential-expression.html
		# spatial RNA
		coords = meta %>%
			select(x,y) %>%
			set_rownames(meta$barcode)
		nUMIs = colSums(input)
		puck = SpatialRNA(coords, input, nUMIs)

		# cell type (label actually? since DE is across all label) specific
		cell_types  = as.factor(meta$cell_type)
		cell_type_threshold = 125 # default
		if (length(unique(cell_types)) == 1){
			cell_types = as.factor(sample(paste0('Celltype', 1:2), length(cell_types), replace=T))
			gene_features$CellTypeB_is_selected = 1
			gene_features$CellTypeA_is_selected = 0
			cell_type_threshold = 0
		}

		names(cell_types) = meta$barcode
		reference = Reference(input, cell_types, nUMIs)
		celltypeA_cols = colnames(gene_features[grepl('CellTypeA', colnames(gene_features)) & grepl('_is_selected', colnames(gene_features))])
		celltypeB_cols = colnames(gene_features[grepl('CellTypeB', colnames(gene_features)) & grepl('_is_selected', colnames(gene_features))])

		if (length(celltypeA_cols) > 1){	
			celltypeA_profiles = rowMeans(gene_features[,celltypeA_cols])
			celltypeB_profiles = rowMeans(gene_features[,celltypeB_cols])	
		} else {
			celltypeA_profiles = gene_features[,celltypeA_cols]
			celltypeB_profiles = gene_features[,celltypeB_cols]
		}
		cell_type_profiles = cbind(celltypeA_profiles, celltypeB_profiles)
		colnames(cell_type_profiles) = c('CellTypeA', 'CellTypeB')
		rownames(cell_type_profiles) = gene_features$Gene

		rctd_res = create.RCTD(puck, reference, cell_type_profiles = cell_type_profiles, max_cores = 2)
		rctd_res = run.RCTD(rctd_res, doublet_mode = 'doublet')
		exp_var = exvar.point.density(rctd_res, rownames(coords), coords[,c('x','y')], radius = 500)
		rctd_res = run.CSIDE.single(rctd_res, exp_var, fdr = 1, log_fc_thresh = 0, gene_threshold = 0, test_genes_sig=F, cell_type_threshold = cell_type_threshold)
		
		# custom run de test
		gene_list_tot = spacexr:::filter_genes(puck, threshold = 0)
		X2 = build.designmatrix.single(rctd_res, exp_var)
		gene_fits = rctd_res@de_results$gene_fits
		cell_types_present = unique(meta$cell_type)
		output_res = data.frame()
		for (cell_type in cell_types_present){
			output_res %<>% rbind(., 
				spacexr:::find_sig_genes_individual(cell_type, cell_types_present, gene_fits,gene_list_tot, X2, params_to_test = 2, fdr = 1, p_thresh = 1,
				log_fc_thresh = 0)$all_genes %>% mutate(cell_type = cell_type, gene = gene_list_tot, p_val_adj = p.adjust(p_val, method = 'BH')))
		}
		output_res %<>%
			select(cell_type, gene, Z_score, p_val, p_val_adj) %>%
			dplyr::rename(stat = Z_score)
	} else if (de_method == 'trendsceek') {
		## only works with spatstat 1.64.1
		## from Giotto https://github.com/RubD/Giotto/blob/HEAD/R/spatial_genes.R
		coords = meta %>%
			select(x,y) %>%
			set_rownames(meta$barcode)

		input = input[rowSums(input) > 0,]
		pp = trendsceek::pos2pp(coords)
		pp = trendsceek::set_marks(pp, as.matrix(input), log.fcn = log10)
		pp[["marks"]] = pp[["marks"]] + 1e-7
		ncores = parallel::detectCores()

		# trendsceektest = trendsceek::trendsceek_test(pp, nrand = 100, ncores=ncores)
		##get rstats
		bp_param = BiocParallel::MulticoreParam(workers = ncores)
		marx = pp[['marks']]
		all_nfeats = ncol(marx)
		all_feats = colnames(marx)

		nfeats = length(selected_features)
		feats = selected_features
		alpha_env = 0.1 / all_nfeats
		alpha_bh = 0.05
		alpha_nom_early = (alpha_bh * 4) / ifelse(all_nfeats >= 500, 10, 1)

		tstat_list = BiocParallel::bplapply(1:nfeats, trendsceek:::calc_trendstats, BPPARAM = bp_param, pp = pp, n.rand = 100, alpha_env = alpha_env, alpha_nom_early = alpha_nom_early)
		names(tstat_list) = feats

		##get supinum stats
		supstats_list = trendsceek:::tstat2supstat(tstat_list)
		trendsceektest = trendsceek:::supstats_list2wide(supstats_list)

		tests = c('Emark', 'Vmark', 'markcorr', 'markvario')
		output_res = data.frame()
		for (test in tests){
			tmp_df = data.frame(
				'gene' = trendsceektest$gene,
				'stat' = trendsceektest[,paste0('max.env.rel.dev_', test)],
				'p_val' = trendsceektest[,paste0('min.pval_', test)],
				'p_val_adj' = trendsceektest[,paste0('p.bh_', test)],
				'test' = test
			)
			output_res %<>% rbind(tmp_df)
		}
	} else if (de_method == 'spark') {
		location = meta %>% select(x, y)
		spark = CreateSPARKObject(counts = input, 
								location = location,
								percentage = 0,
								min_total_counts = 0)
		spark@lib_size = apply(spark@counts, 2, sum)
		ncores = parallel::detectCores()
		spark = spark.vc(spark, 
						covariates = NULL, 
						lib_size = spark@lib_size, 
						num_core = ncores, 
						verbose = F)
		spark = spark.test(spark, 
						check_positive = T,
						verbose = F)
		output_res = data.frame(
			'gene' = rownames(spark@res_mtest),
			'p_val' = spark@res_mtest$combined_pvalue
		)
		output_res$p_val_adj = p.adjust(output_res$p_val, method = 'BH')
	} else if (de_method == 'sparkx') {
		location = meta %>% select(x, y)
		sparkX = sparkx(input, location, numCores=1, option='mixture')
		output_res = data.frame(
			'gene' = rownames(sparkX$res_mtest),
			'p_val' = sparkX$res_mtest$combinedPval
		)
		output_res$p_val_adj = p.adjust(output_res$p_val, method = 'BH')
	} else if (de_method == 'heartsvg') {
		# follow default params
		scale = T
		qh = 0.985
		noise = T
		padj_m = 'holm'

		# custom create data.frame with first 2 columns as x,y and other cols as genes
		data = t(input)
		coords = meta %>% 
			select(x,y) %>% set_rownames(meta$barcode) %>% set_colnames(c('row', 'col'))
		coords = coords[rownames(data),]
		data = cbind(coords, data)

		scale.count=function(data,qh=0.985){
			counts_mat=apply(data, 1, function(y){z=y/max(y,na.rm = T)})
			counts_mat=apply(counts_mat, 2, function(y)
				{y[which(y>quantile(y,qh,na.rm=T))]=quantile(y,qh,na.rm=T);
					y[which(y<quantile(y,0.25,na.rm=T))]=0;y})
			t(counts_mat)
		}

		if (scale == T) {
			new = cbind(data[, 1:2], scale.count(data[, -c(1:2)], qh))
		}
		if (scale == F) {
			new = data
		}
		new = as.data.frame(new)
		locus_in = new[, c(1, 2)]
		counts = new[, -c(1, 2)]
		l1 = "row"
		l2 = "col"
		z1_group = ceiling(log(diff(range(locus_in[l1]))))
		z2_group = ceiling(log(diff(range(locus_in[l2]))))

		locus_in1 = data.frame(locus_in, n.row = cut(x = as.matrix(locus_in[l1]), breaks = seq(from = min(locus_in[l1]), to = max(locus_in[l1]), length = z1_group + 1), labels = 1:z1_group, include.lowest = T, right = T))
		locus_in2 = data.frame(locus_in, n.col = cut(x = as.matrix(locus_in[l2]), breaks = seq(from = min(locus_in[l2]), to = max(locus_in[l2]), length = z2_group + 1), labels = 1:z2_group, include.lowest = T, right = T))
		new1 = cbind(locus_in1[colnames(locus_in1) != l1], counts)
		new2 = cbind(locus_in2[colnames(locus_in2) != l2], counts)
		colnames(new1)[1:2] = c("x", "coor.z")
		colnames(new2)[1:2] = c("x", "coor.z")
		new_row = aggregate(new1[, -c(1:2)], list(new1$x, new1$coor.z),
							mean)
		new_col = aggregate(new2[, -c(1:2)], list(new2$x, new2$coor.z),
							mean)
		zero.p = apply(new[, -c(1:2)], 2, function(y) {
			sum(y != 0, na.rm = T)/length(y)
		})
		mean = apply(new[, -c(1:2)], 2, function(y) {
			mean(y[y != 0], na.rm = T)
		})
		sum = data.frame(gene = names(mean), zero.p, mean)
		row_t = apply(new_row[, -c(1:2)], 2, function(y) {
			z = ifelse(sum(y != 0) == 0, 1, Box.test(y, lag = z1_group)$p.value)
			z
		})
		col_t = apply(new_col[, -c(1:2)], 2, function(y) {
			z = ifelse(sum(y != 0) == 0, 1, Box.test(y, lag = z2_group)$p.value)
			z
		})
		new_x = aggregate(new[, -c(1:2)], list(new$row), mean)
		new_y = aggregate(new[, -c(1:2)], list(new$col), mean)
		x_t = apply(new_x[, -1], 2, function(y) {
			Box.test(y, lag = z1_group)$p.value
		})
		y_t = apply(new_y[, -1], 2, function(y) {
			Box.test(y, lag = z2_group)$p.value
		})
		test = data.frame(row_t, col_t, x_t, y_t, gene = names(row_t))
		test[, 1:(ncol(test) - 1)] = sapply(test[, 1:(ncol(test) - 1)], function(y) {
			z = ifelse(is.na(y) == T, 0.999, y)
		})
		test$min = apply(test[, 1:(ncol(test) - 1)], 1, function(y) {
			poolr::stouffer(y)$p
		})
		test$adj_min = p.adjust(test$min, method = padj_m)
		test = merge(test, sum, by = "gene")
		test = subset(test, zero.p > 0)
		test$zero.p = test$zero.p/max(test$zero.p)
		test$mean = test$mean/max(test$mean)
		test$c2 = (test$zero.p + test$mean)/2
		if (noise == T) {
			data.table::setorder(test, adj_min, min, -zero.p, -mean)
			a1 = c("gene", "min", "adj_min")
			test2 = test[a1]
			b1 = c("gene", "pval", "p_adj")
			colnames(test2) = b1
			test2$rank = 1:nrow(test2)
		}
		if (noise == F) {
			data.table::setorder(test, adj_min, min)
			a1 = c("gene", "min", "adj_min")
			test2 = test[a1]
			b1 = c("gene", "pval", "p_adj")
			colnames(test2) = b1
			test2$rank = 1:nrow(test2)
		}

		test2$p_val_adj = p.adjust(test2$pval, method = 'BH')
		output_res = test2 %>%
			mutate(
				p_val = pval
			) %>%
			select(gene, p_val, p_val_adj) %>%
			set_rownames(NULL)
	} else if (de_method == 'spade') {
		labels = unique(meta$label)
		
		cells1 = meta %>% filter(label == labels[1]) %>% pull(barcode)
		counts1 = input[, cells1]
		info1 = meta %>% filter(label == labels[1]) %>% select(x,y)

		cells2 = meta %>% filter(label == labels[2]) %>% pull(barcode)
		counts2 = input[, cells2]
		info2 = meta %>% filter(label == labels[2]) %>% select(x,y)
		
		stable1 = stabilize(as.matrix(counts1))
		info1$total_counts = colSums(counts1)
		info_py1 = r_to_py(info1)
		df_py1 = r_to_py(as.data.frame(stable1))
		regressed1 = regress_out(info_py1, df_py1, "np.log(total_counts)")
		rownames(regressed1) = rownames(counts1)
		colnames(regressed1) = colnames(counts1)
		regressed1 = as.matrix(regressed1)

		stable2 = stabilize(as.matrix(counts2))
		info2$total_counts = colSums(counts2)
		info_py2 = r_to_py(info2)
		df_py2 = r_to_py(as.data.frame(stable2))
		regressed2 = regress_out(info_py2, df_py2, "np.log(total_counts)")
		rownames(regressed2) = rownames(counts2)
		colnames(regressed2) = colnames(counts2)
		regressed2 = as.matrix(regressed2)

		# maunal run output_res = SPADE_DE(regressed1, regressed2, info1, info2) with parallel lol
		mode="Shape&Strength"
		ED1 = as.matrix(dist(info1))
		lrang1 = ComputeGaussianPL(ED1, compute_distance=FALSE)
		ED2 = as.matrix(dist(info2))
		lrang2 = ComputeGaussianPL(ED2, compute_distance=FALSE)
		ncores = parallel::detectCores()
		if (ncores > 1){
			lapply_fun = parallel::mclapply
		} else {
			lapply_fun = lapply
		}
		genes = rownames(regressed1)
		para = data.frame(do.call(rbind, lapply_fun(genes, function(gene){
			message(gene)
			y1 = regressed1[gene,]
			y2 = regressed2[gene,]

			re1 = optimize(lengthscale_fit, c(lrang1[3],lrang1[8]), location=info1, y=y1, tol=1)
			re2 = optimize(lengthscale_fit, c(lrang2[3],lrang2[8]), location=info2, y=y2, tol=1)

			Est1 = Delta_fit(location=info1, y=y1, L=re1$minimum)
			Est2 = Delta_fit(location=info2, y=y2, L=re2$minimum)

			delta1 = Est1$delta
			delta2 = Est2$delta

			Tao_hat1 = Est1$Tao_hat
			Tao_hat2 = Est2$Tao_hat

			if (mode=="Shape&Strength"){
				LL10 = LL_DE(delta=delta2, location=info1, L=re2$minimum, y=y1)
				LL20 = LL_DE(delta=delta1, location=info2, L=re1$minimum, y=y2)
			}
			if (mode=="Shape"){
				LL10 = LL_DE(delta=delta1, location=info1, L=re2$minimum, y=y1)
				LL20 = LL_DE(delta=delta2, location=info2, L=re1$minimum, y=y2)
			}
			if (mode=="Strength"){
				LL10 = LL_DE(delta=delta2, location=info1, L=re1$minimum, y=y1)
				LL20 = LL_DE(delta=delta1, location=info2, L=re2$minimum, y=y2)
			}
			data.frame(
				gene = gene,
				theta_Gau1 = re1$minimum,
				theta_Gau2 = re2$minimum,
				Gamma1 = Tao_hat1,
				Gamma2 = Tao_hat2,
				logLik11 = -re1$objective,
				logLik21 = -re2$objective,
				logLik10 = LL10,
				logLik20 = LL20
			)
		})))

		para$Diff = 2*(para$logLik11 + para$logLik21 - para$logLik10 - para$logLik20)
		para$p_val = pchisq((para$Diff),df=1,lower.tail=FALSE)
		para$p_val_adj = p.adjust(para$p_val, method = "BH")
		output_res = para %>%
			dplyr::select(gene, p_val, p_val_adj)
	} else if (de_method %in% c('scran', 'mast', 'binSpect_kmeans', 'binSpect_rank')){
		coords_df = meta %>%
			select(x,y) %>%
			set_rownames(meta$barcode)
		giotto_instructions = createGiottoInstructions(python_path = giotto_python_path)
		input = input[rowSums(input) > 0,]
		giotto_obj = createGiottoObject(input, cell_metadata = meta, spatial_locs = coords_df, instructions = giotto_instructions)
		giotto_obj %<>% normalizeGiotto()

		if (de_method == 'scran'){
			gene_markers = findScranMarkers_one_vs_all(
				gobject = giotto_obj,
				expression_values = 'normalized',
				cluster_column = 'label',
				pval = 1,
				logFC = 0
			)
			output_res = gene_markers %>%
				as.data.frame() %>%
				select(genes, p.value) %>%
				set_colnames(c('gene', 'p_val')) 
			output_res$p_val_adj = p.adjust(output_res$p_val, method='BH')	
		} else if (de_method == 'mast') {
			gene_markers = findMastMarkers(
				gobject = giotto_obj,
				expression_values = 'normalized',
				cluster_column = 'label',
				group_1 = 'Perturbation_1',
				group_1_name = 'Perturbation_1',
				group_2 = 'Perturbation_2',
				group_2_name = 'Perturbation_2'
			)
			output_res = gene_markers %>%
				as.data.frame() %>%
				select(genes, `Pr(>Chisq)`) %>%
				set_colnames(c('gene', 'p_val'))
			output_res$p_val_adj = p.adjust(output_res$p_val, method='BH')
		} else if (de_method == 'binSpect_kmeans'){
			giotto_obj %<>% createSpatialDelaunayNetwork()
			km_spatialgenes = binSpect(giotto_obj, bin_method = 'kmeans')
			output_res = km_spatialgenes %>%
				select(genes, p.value) %>%
				set_colnames(c('gene', 'p_val'))
			output_res$p_val_adj = p.adjust(output_res$p_val, method='BH')	
		} else if (de_method == 'binSpect_rank'){
			giotto_obj %<>% createSpatialDelaunayNetwork()
			km_spatialgenes = binSpect(giotto_obj, bin_method = 'rank')
			output_res = km_spatialgenes %>%
				select(genes, p.value) %>%
				set_colnames(c('gene', 'p_val'))
			output_res$p_val_adj = p.adjust(output_res$p_val, method='BH')	
		} 
	}
	return (output_res)
}

