## import dependencies 
import argparse
import pandas as pd
import os
import numpy as np
import scanpy as sc
import math
import sys
import SpaGFT as spg

# set working directory
git_dir = os.path.expanduser("~/git/vespucci")
python_dir = git_dir + "/python"
os.chdir(python_dir)
sys.path.append(python_dir)

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

# from https://spagft.readthedocs.io/en/latest/spatial/lymphnode_tutorial.html#4.-Function:-identify-spatially-variable-genes
adata = sc.read(args.anndata_input_filename)
adata.obs.rename(columns = {'x': 'array_col', 'y':'array_row'}, inplace = True)
adata.obsm['spatial'] = adata.obs[['array_col', 'array_row']]
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
# determine the number of low-frequency FMs and high-frequency FMs
(ratio_low, ratio_high) = spg.gft.determine_frequency_ratio(adata,ratio_neighbors=1)
gene_df = spg.detect_svg(adata, spatial_info=['array_row', 'array_col'], ratio_low_freq=ratio_low, ratio_high_freq=ratio_high, ratio_neighbors=1,filter_peaks=True, S=6)
output_res = gene_df[['gft_score','pvalue']].drop_duplicates()
output_res['gene'] = output_res.index
output_res = output_res.reset_index(drop=True)
output_res.rename(columns = {'pvalue': 'p_val'}, inplace = True)
output_res.to_csv(args.output_filename)