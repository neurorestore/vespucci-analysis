setwd("~/git/vespucci-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-align-spaceranger.R')
grid = read.delim("sh/grids/regen/align-spaceranger.txt")
for (param_name in colnames(grid))
parser$add_argument(paste0('--', param_name),
type = typeof(grid[[param_name]]))

args = parser$parse_args()
print(args)

# detect system
set.seed(1000)

library(tidyverse)
library(magrittr)

output_dir = args$output_dir
if (!dir.exists(output_dir)) {
	dir.create(output_dir, recursive = TRUE)
}

setwd(output_dir)
spaceranger_path = 'spaceranger-2.1.1/spaceranger'
fun_call = paste0(spaceranger_path, ' count')
input = paste0('--id=', args$id)
output = paste0('--output_dir=', args$output_dir)
transcriptome = paste0('--transcriptome=', args$ref_genome)
fastqs = paste0('--fastqs=', args$fastq_dir)
sample = paste0('--sample=', args$id)
image = paste0('--image=', args$image_file)
slide = paste0('--slide=', args$slide_ref)
area = paste0('--area=', args$area)
other_args = paste0('--reorient-images true')

full_call = paste(
	fun_call, 
	input, 
	transcriptome,
	fastqs,
	sample,
	image,
	slide,
	area,
	other_args
)

# Run the call
system(full_call)
