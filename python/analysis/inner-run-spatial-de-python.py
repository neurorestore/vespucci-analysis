## import dependencies 
import argparse
import pandas as pd
import os
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import tensorflow as tf
import sys
# import SpatialDE
import NaiveDE
import squidpy as sq
# import smashpy
# import shap
# from somde import SomNode

# set working directory
git_dir = os.path.expanduser("~/git/vespucci")
python_dir = git_dir + "/python"
os.chdir(python_dir)
sys.path.append(python_dir)

from analysis.spatialde_modules import *
from analysis.smashpy_modules import *

### dynamically build CLI
parser = argparse.ArgumentParser()
## build the CLI
grid_file = git_dir + '/sh/grids/simulations/calculate-spatial-de-python.txt'
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
	output_res = pd.DataFrame(data={'genes': genes, "p_val":pvals, "p_val_adj":pvals_adj})
	output_res.to_csv(args.output_filename)

elif args.de_method == 'spatialDE':
	adata = sc.read(args.anndata_input_filename)
	counts = pd.DataFrame(adata.raw.X.toarray(), columns=adata.var['Gene'], index = adata.obs['Cell'])
	meta = pd.DataFrame({'x': adata.obs['x'], 'y': adata.obs['y'], 'total_counts': adata.obs['nCount_originalexp']}, index = adata.obs['Cell'])

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
	adata.obsm['spatial'] = pd.DataFrame({'x': adata.obs['x'], 'y': adata.obs['y']}, index = adata.obs['Cell'])
	sq.gr.spatial_neighbors(adata)
	genes = adata.var['Gene']
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