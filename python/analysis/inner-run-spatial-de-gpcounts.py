## import dependencies 
import argparse
import pandas as pd
import os
import numpy as np
import scanpy as sc
import math
import gpflow ## only work with 2.0.5
import tensorflow as tf
from GPcounts.GPcounts_Module import Fit_GPcounts
import sys

# set working directory
git_dir = os.path.expanduser("~/git/vespucci-analysis/")
python_dir = git_dir + "/python"
os.chdir(python_dir)
sys.path.append(python_dir)

### dynamically build CLI
parser = argparse.ArgumentParser()
## build the CLI
grid_file = git_dir + '/sh/grids/simulations/calculate-spatial-de-gpcounts.txt'
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

adata = sc.read(args.anndata_input_filename)
data = pd.DataFrame(np.transpose(adata.raw.X.toarray()), columns=adata.obs['Cell'], index=adata.var['Gene'])
spatial_locations = pd.DataFrame({'x': adata.obs['x'], 'y': adata.obs['y']}, index = adata.obs['Cell'])
likelihood = 'Negative_binomial'
gp_counts = Fit_GPcounts(spatial_locations, data)
log_likelihood_ratio = gp_counts.One_sample_test(likelihood)
output_res = gp_counts.calculate_FDR(log_likelihood_ratio)
output_res.rename(columns = {'p_value':'p_val', 'log_likelihood_ratio':'stat', 'q_value': 'p_val_adj'}, inplace = True)
output_res['gene'] = output_res.index
output_res = output_res[['gene', 'stat', 'p_val', 'p_val_adj']]
output_res.to_csv(args.output_filename)