setwd("~/git/vespucci-analysis")
options(stringsAsFactors = FALSE)
library(argparse)
library(Seurat)
library(tidyverse)
library(magrittr)
library(Vespucci)

# dynamically parse arguments

parser = ArgumentParser(prog = 'inner-run-vespucci.R')
grid = read.delim("sh/grids/simulations/run-vespucci-grid.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name), type = typeof(grid[[param_name]]))

args = parser$parse_args()
print(args)

sc = readRDS(args$input_filename)
expr = GetAssayData(sc, slot='counts')
meta = sc@meta.data

ves_res = run_vespucci(
    input = expr,
    meta = meta,
    max_barcodes = 100
)

saveRDS(ves_res, args$ves_output_filename)

auc_res = ves_res$spatial_auc_result$aucs
de_res = ves_res$de_feature_result

gene_features = sc@assays$originalexp@meta.features
min_pval = min(de_res$p_val[de_res$p_val > 0])

perturb_1_cols = colnames(gene_features)[grepl('Perturbation_1', colnames(gene_features)) & grepl('_is_selected', colnames(gene_features))]
perturb_2_cols = colnames(gene_features)[grepl('Perturbation_2', colnames(gene_features)) & grepl('_is_selected', colnames(gene_features))]
gene_features$de_effect = rowSums(gene_features[,c(perturb_1_cols, perturb_2_cols)]) > 0
gene_features$gene = gene_features$Gene

output = gene_features %>%
    left_join(de_res, by = 'gene') %>%
    filter(!is.na(p_val)) %>%
    dplyr::select(gene, de_effect, logFC, p_val, p_val_adj)

output %<>% 
    mutate(
        auc_de = ifelse(p_val > 0, -log10(p_val), -log10(min_pval)),
        auc_de_binary = p_val_adj < 0.05,
        truth = as.numeric(de_effect)
    )

output %>%
    dplyr::select(auc_de_binary, truth) %>%
    table()

auroc_val = 
    tryCatch({
        roc_res = roc(
            predictor = output$auc_de,
            response = factor(output$truth, levels = c(0,1))
        )
        auc(roc_res)[1]
    },
    error = function(e){
        message(e)
        return (-1)
    })

output %<>%
    mutate(
        auc_de_binary_int = as.integer(auc_de_binary)
    )

auprc_vals = 
    tryCatch({
        pr_res = pr.curve(scores.class0 = output$auc_de[!is.na(output$auc_de) & output$truth==1], scores.class1 = output$auc_de[!is.na(output$auc_de) & output$truth==0])
        list(
            'auprc_integral' = pr_res$auc.integral,
            'auprc_davis_goadrich' = pr_res$auc.davis.goadrich
        )
    },
    error = function(e){
        message(e)
        return(
            list(
                'auprc_integral' = -1,
                'auprc_davis_goadrich' = -1
            )
        )
    })

tp_genes = output %>%
    filter(auc_de_binary_int == truth, truth == 1) %>%
    pull(gene)
tn_genes = output %>%
    filter(auc_de_binary_int == truth, truth == 0) %>%
    pull(gene)
fp_genes = output %>%
    filter(auc_de_binary_int != truth, auc_de_binary == 1) %>%
    pull(gene)
fn_genes = output %>%
    filter(auc_de_binary_int != truth, auc_de_binary != 1) %>%
    pull(gene)

tp_size = length(tp_genes)
tn_size = length(tn_genes)
fp_size = length(fp_genes)
fn_size = length(fn_genes)

sensitivity = tp_size/(tp_size + fn_size)
specificity = tn_size/(tn_size + fp_size)
ppv = tp_size/(tp_size + fp_size)
npv = tn_size/(tn_size + fn_size)

acc = (tp_size + tn_size)/(nrow(output))

mcc = mcc(output$auc_de_binary_int, output$truth)

stats = data.frame(
    auroc = auroc_val,
    auprc_integral = auprc_vals[['auprc_integral']],
    auprc_davis_goadrich = auprc_vals[['auprc_davis_goadrich']],
    tp = tp_size,
    tn = tn_size,
    fp = fp_size,
    fn = fn_size,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    acc = acc,
    mcc = mcc
)

stats[is.na(stats)] = 0

output_list = list(
    'nebula_res' = output,
    'stats' = stats
)
saveRDS(output_list, args$de_output_filename)
