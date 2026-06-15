## import dependencies 
import argparse
import pandas as pd
import os
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
from scipy import optimize
import random, torch
import tensorflow as tf
import sys
import SpatialDE
import NaiveDE
import squidpy as sq
from somde import SomNode, dyn_de, get_mll_results, stabilize, regress_out
# from scGCO import normalize_count_cellranger, create_graph_with_weight
from scGCO import *
# import smashpy
# import shap
# from somde import SomNode

# set working directory
git_dir = os.path.expanduser("~/git/vespucci-analysis/")
python_dir = git_dir + "/python"
os.chdir(python_dir)
sys.path.append(python_dir)

from analysis.spatialde_modules import *
from analysis.spanve_modules import *
# from analysis.smashpy_modules import *
# import time
# import resource

### dynamically build CLI
parser = argparse.ArgumentParser()
## build the CLI
grid_file = git_dir + '/sh/grids/real_data/calculate-spatial-de-python.txt'
grid = pd.read_csv(grid_file, sep='\t')
for arg_name in list(grid):
	param_name = '--' + arg_name
	param_dtype = str(grid[arg_name].dtype)
	# convert to pandas
	param_type = {'object': str,
				  'int64': int,
				  'float64': float,
				  'bool': bool
				  }[param_dtype]
	parser.add_argument(param_name, type=param_type)

# parse all arguments
args = parser.parse_args()
print(args)

# args.anndata_input_filename = '/work/upcourtine/vespucci/simulations/debug/python_input/anndata/test_sample.h5ad'

time_start = time.perf_counter()
if args.de_method == 'spacgn':
	adata=sc.read(args.anndata_input_filename)
	adata.var_names_make_unique()

	# what should we separate this by? label or cell type or all as one?
	# cell type because SVG only works within each domain (pred)
	adata.obs["pred"]=adata.obs['label'].astype('category')
	adata.obs["x_pixel"]=adata.obs["x"]
	adata.obs["y_pixel"]=adata.obs["y"]
	x_pixel=adata.obs["x"].tolist()
	y_pixel=adata.obs["y"].tolist()
	#Convert sparse matrix to non-sparse
	adata.X=adata.raw.X.toarray()
	#adata.adata=adata
	sc.pp.log1p(adata)

	#Use domain 0 as an example
	target=adata.obs['pred'].unique()[0]
	neigbour=adata.obs['pred'].unique()[1]
	#Set filtering criterials
	min_in_group_fraction=0.8
	min_in_out_group_ratio=1
	min_fold_change=1.5

	adata.obs["target"]=((adata.obs['pred']==target)*1).astype('category')
	sc.tl.rank_genes_groups(adata, groupby="target",reference="rest", n_genes=adata.shape[1],method='wilcoxon')
	pvals=[i[0] for i in adata.uns['rank_genes_groups']["pvals"]]
	pvals_adj=[i[0] for i in adata.uns['rank_genes_groups']["pvals_adj"]]
	genes=[i[1] for i in adata.uns['rank_genes_groups']["names"]]
	if args.datatype == 'real':
		genes =adata.var.features[[int(x) for x in genes]]
	output_res = pd.DataFrame(data={'genes': genes, "p_val":pvals, "p_val_adj":pvals_adj})
	output_res.to_csv(args.output_filename)

elif args.de_method == 'spatialDE':
	adata = sc.read(args.anndata_input_filename)
	if args.datatype == 'sim':
		counts = pd.DataFrame(adata.raw.X.toarray(), columns=adata.var['Gene'], index = adata.obs['Cell'])
		meta = pd.DataFrame({'x': adata.obs['x'], 'y': adata.obs['y'], 'total_counts': adata.obs['nCount_originalexp']}, index = adata.obs['Cell'])
	if args.datatype == 'real':
		counts = pd.DataFrame(adata.raw.X.toarray(), columns=adata.var['features'], index = adata.obs['barcode'])
		meta = pd.DataFrame({'x': adata.obs['x'], 'y': adata.obs['y'], 'total_counts': adata.obs['nCount_RNA']}, index = adata.obs['barcode'])

	# meta = meta.sample(1000)
	counts = counts.loc[meta.index]
	norm_expr = NaiveDE.stabilize(counts.T).T
	resid_expr = NaiveDE.regress_out(meta, norm_expr.T, 'np.log(total_counts)').T
	# resid_expr = resid_expr.sample(n=100, axis=1, random_state=1)
	X = meta[['x', 'y']]
	l_min, l_max = get_l_limits(X)
	kernel_space = {
		'SE': np.logspace(np.log10(l_min), np.log10(l_max), 10),
		'const': 0
	}
	results = []
	result = const_fits(resid_expr)
	result['l'] = np.nan
	result['M'] = 2
	result['model'] = 'const'
	results.append(result)

	US_mats = []
	# t0 = time()
	for ii, lengthscale in enumerate(kernel_space['SE']):
		print(ii)
		K = SE_kernel(X, lengthscale)
		U, S = factor(K)
		gower = gower_scaling_factor(K)
		UT1 = get_UT1(U)
		US_mats.append({
			'model': 'SE',
			'M': 4,
			'l': lengthscale,
			'U': U,
			'S': S,
			'UT1': UT1,
			'Gower': gower
		})

	n_models = len(US_mats)

	for i, cov in enumerate(tqdm(US_mats, desc='Models: ')):
		result = lengthscale_fits(resid_expr, cov['U'], cov['UT1'], cov['S'], cov['Gower'])
		result['l'] = cov['l']
		result['M'] = cov['M']
		result['model'] = cov['model']
		results.append(result)

	n_genes = resid_expr.shape[1]
	results = pd.concat(results, sort=True).reset_index(drop=True)
	results['BIC'] = -2 * results['max_ll'] + results['M'] * np.log(results['n'])

	mll_results = get_mll_results(results)
	# Perform significance test
	mll_results['pval'] = 1 - stats.chi2.cdf(mll_results['LLR'], df=1)

	output_res = mll_results[['g', 'LLR', 'pval']]
	output_res.rename(columns = {'g':'gene', 'pval': 'p_val'}, inplace = True)
	output_res.to_csv(args.output_filename)

elif args.de_method == 'squidpy':
	adata = sc.read(args.anndata_input_filename)
	adata.X=adata.raw.X.toarray()
	if args.datatype == 'sim':
		adata.obsm['spatial'] = pd.DataFrame({'x': adata.obs['x'], 'y': adata.obs['y']}, index = adata.obs['Cell'])
		genes = adata.var['Gene']
	elif args.datatype == 'real':
		adata.obsm['spatial'] = pd.DataFrame({'x': adata.obs['x'], 'y': adata.obs['y']}, index = adata.obs['barcode'])
		genes = adata.var['features']
	sq.gr.spatial_neighbors(adata)
	sq.gr.spatial_autocorr(
		adata,
		mode="moran",
		genes=genes,
		n_perms=100,
		n_jobs=1,
	)
	tmp_res = adata.uns["moranI"]
	output_res = pd.concat([
		pd.DataFrame(data={'genes': tmp_res.index, "p_val":tmp_res['pval_norm'], "p_val_adj":tmp_res['pval_norm_fdr_bh'], 'type':'normality'}),
		pd.DataFrame(data={'genes': tmp_res.index, "p_val":tmp_res['pval_sim'], "p_val_adj":tmp_res['pval_sim_fdr_bh'], 'type':'permutation'}),
		pd.DataFrame(data={'genes': tmp_res.index, "p_val":tmp_res['pval_z_sim'], "p_val_adj":tmp_res['pval_z_sim_fdr_bh'], 'type':'normal_approx_permutation'}),
		], axis=0)
	output_res.to_csv(args.output_filename)

elif args.de_method == 'smash':
	sm = smashpy.smashpy()
	adata = sc.read(args.anndata_input_filename)
	adata.X = adata.raw.X.toarray()
	adata.obs['annotation'] = ('label' + adata.obs['label'].astype('str')).astype('category')
	sm.data_preparation(adata)
	sm.DNN(adata, group_by="annotation", model=None, balance=True, verbose=True, save=False)

	# manual run sm.run_shap(adata, group_by="annotation", model=None, verbose=True, pct=0.001, restrict_top=("local", 20))

	X = np.array(adata.X)
	y = np.array(adata.obs['annotation'].tolist())

	myDict = {}
	for idx, c in enumerate(adata.obs['annotation'].cat.categories):
		myDict[c] = idx

	labels = []
	for l in adata.obs['annotation'].tolist():
		labels.append(myDict[l])

	labels = np.array(labels)
	y = labels
	SEED = 42
	pct = 0.05
	num_classes = len(np.unique(adata.obs['annotation']))
	model = loadDNNmodel(adata, num_classes)
	model.compile(loss=losses.categorical_crossentropy,
				  optimizer=optimizers.Adam(learning_rate=0.001, amsgrad=False),
				  metrics=['accuracy', 'AUC', 'Precision', 'Recall'])
	model.load_weights('weights/best_model_%s.h5'%'annotation')
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=SEED, stratify=y)
	X_tr1, x_val1, Y_tr1, y_val1 = train_test_split(X_train, y_train, test_size=(1-pct), random_state=SEED, stratify=y_train)
	X_tr2, x_val2, Y_tr2, y_val2 = train_test_split(X_test, y_test, test_size=pct, random_state=SEED, stratify=y_test)
	explainer = shap.DeepExplainer(model, X_tr1)
	shap.explainers._deep.deep_tf.op_handlers["AddV2"] = shap.explainers._deep.deep_tf.passthrough
	shap_values = explainer.shap_values(x_val2)
	output_res = pd.concat([
		pd.DataFrame({'gene': adata.var['Gene'], 'label': 'label1', 'shap_value': shap_values[0][0,]}),
		pd.DataFrame({'gene': adata.var['Gene'], 'label': 'label2', 'shap_value': shap_values[0][1,]})],
		axis=0
		)
	output_res.to_csv(args.output_filename)

elif args.de_method == 'spatialDE2':
	adata = sc.read(args.anndata_input_filename)
	# adata = adata[:,:10]
	adata.obsm['spatial'] = adata.obs[['x', 'y']]
	svg_full, _ = SpatialDE.test(adata, omnibus=True)
 
	output_res = svg_full[['gene', 'kappa', 'pval']]
	output_res.rename(columns = {'pval': 'p_val'}, inplace = True)
	output_res.to_csv(args.output_filename)

elif args.de_method == 'somde':
	# from https://pypi.org/project/somde/
	adata = sc.read(args.anndata_input_filename)
	X = adata.obs[['x', 'y']].values.astype(np.float32)
	som = SomNode(X, 20)
	if args.datatype == 'real': # for real data
		df = pd.DataFrame.sparse.from_spmatrix(adata.X, columns=adata.var.features, index=adata.obs.index).transpose()
	elif args.datatype == 'sim': # for simulations
		df = pd.DataFrame(adata.X, columns=adata.var.Gene, index=adata.obs.index).transpose()
	ndf,ninfo = som.mtx(df)
	nres = som.norm()
 	# expression_matrix = ndf
    # phi_hat, _ = optimize.curve_fit(lambda mu, phi: mu + phi * mu ** 2, expression_matrix.mean(1), expression_matrix.var(1))
	# np.log(expression_matrix + 1. / (2 * phi_hat[0]))
	# som.run()
	# from https://github.com/WhirlFirst/somde/blob/1f015b4ee90100fadfae04631e8855a0fac2c9bf/somde-python/som.py since som.run() throws error
	X1 = som.ninfo[['x','y']].values.astype(float)
	l_min, l_max = get_l_limits(X1)
	kernel_space = {
		'SE': np.logspace(np.log10(l_min), np.log10(l_max), 10),
		'const': 0
	}
	results = dyn_de(X1, som.nres, kernel_space)
	mll_results = get_mll_results(results)
	# Perform significance test
	mll_results['pval'] = 1 - stats.chi2.cdf(mll_results['LLR'], df=1)
	output_res = mll_results[['g', 'LLR', 'pval']].drop_duplicates().reset_index(drop=True)
	output_res.rename(columns = {'g': 'gene', 'pval': 'p_val'}, inplace = True)
	output_res.to_csv(args.output_filename)

elif args.de_method == 'spanve':
	# from https://github.com/zjupgx/Spanve/tree/main
	adata = sc.read(args.anndata_input_filename)
	adata.obsm['spatial'] = adata.obs[['x', 'y']]
	spanve = Spanve(adata)
	spanve.fit()
	output_res = pd.DataFrame(spanve.result_df[['ent', 'pvals']]).reset_index(drop=True)
	output_res.rename(columns = {'pvals': 'p_val'}, inplace = True)
	output_res['gene'] = spanve.result_df.index.values
	output_res = output_res[['gene', 'ent', 'p_val']]
	output_res.to_csv(args.output_filename)
 
elif args.de_method == 'scgco':
	adata = sc.read(args.anndata_input_filename)
	# The required matrix format is the ST data format, a matrix of counts where spot coordinates are row names and the gene names are column names. This default matrix format (.TSV ) is split by tab.
	if args.datatype == 'real':
		data = pd.DataFrame.sparse.from_spmatrix(adata.X, columns=adata.var.features, index=adata.obs.index).sparse.to_dense()
	elif args.datatype == 'sim':
		data = pd.DataFrame(adata.X, columns=adata.var.Gene, index=adata.obs.index)
	# data = data.iloc[1:1000,1:100]
	data = data.loc[:,data.sum(axis=0)>0]
	# data_norm = normalize_count_cellranger(data)
	# normalize_count_cellranger doesn't work so from https://github.com/WangPeng-Lab/scGCO/blob/master/code/scGCO_code/scGCO_source/Preprocessing.py
	normalizing_factor = np.sum(data, axis = 1)/np.median(np.sum(data, axis = 1))
	data_norm = data.divide(normalizing_factor, axis=0)
	data_norm = np.log1p(data_norm)
	locs = np.array(adata.obs[['x', 'y']].iloc[1:1000,:].reset_index(drop=True))
	exp = data_norm.iloc[:,0]
	# needs to replace envs/de-methods-python/.../site-packages/scGCO/Preprocessing.py
	cellGraph = create_graph_with_weight(locs, exp)
	gmmDict = multiGMM(data_norm)
	# replace line "sparse_matrix = nx.to_scipy_sparse_matrix(graph)" to "sparse_matrix = nx.nx.to_scipy_sparse_array(graph)" in envs/de-methods-python/.../site-packages/pysal/lib/weights/weights.py
	result_df = identify_spatial_genes(locs, data_norm, cellGraph, gmmDict)
	best_p_values=[min(i) for i in result_df['p_value']]
	output_res = pd.DataFrame({'gene':result_df.index,'p_val':best_p_values})
	output_res.to_csv(args.output_filename)