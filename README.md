# README

Scripts used to reproduce all figures in _"Identification of perturbation-responsive regions and genes in comparative spatial transcriptomics atlases"_.

---

## Repository layout

```
R/figures/          # one script per extended figure panel (see table below)
R/functions/        # shared utility functions
R/theme.R           # shared ggplot theme
data/               # symlink to the data directory (see Prerequisites)
```

All scripts call `setwd('~/git/vespucci-analysis')` at the top and load data via the `data/` symlink.

---

## Prerequisites

### 1. Clone this repository

```bash
git clone https://github.com/CSOgroup/vespucci-analysis.git  ~/git/vespucci-analysis
```

### 2. Data

Download the data archive from Zenodo (DOI: _to be added_) and extract it to a local directory, then symlink it as `data/` inside this repo:

```bash
ln -s /path/to/extracted/data  ~/git/vespucci-analysis/data
```

The expected directory structure inside `data/` is:

```
data/
├── simulations/
│   ├── objects/                   # Seurat simulation objects (.rds)
│   ├── objects_meta/              # gene + cell metadata per simulation
│   ├── vespucci/                  # Vespucci results on simulations
│   ├── registration/              # registration robustness experiments
│   │   ├── meta/
│   │   └── vespucci/
│   ├── summaries/
│   │   ├── spatial_acc.rds
│   │   ├── simulations-auc-summary.rds
│   │   ├── simulations-pvals-summary.rds
│   │   ├── magellan/
│   │   ├── spatial_res/
│   │   ├── timeit/
│   │   └── de_results/
│   │       ├── all_vespucci_de_auroc_summary.rds
│   │       ├── other_methods_stats.rds
│   │       ├── vespucci_de_auroc_false_summary.rds
│   │       ├── other_methods_false_stats.rds
│   │       └── splatter_de_false_results.rds
│   └── DE_summaries/
│       ├── others/
│       └── others_python/
├── real_data/
│   ├── seurat/                    # processed Seurat objects per dataset
│   │   ├── regen_final.rds        # all-conditions regen Seurat object (main Fig 4)
│   │   └── ...
│   ├── seurat_GO/DE/              # GO-module Seurat objects
│   ├── meta/                      # coordinate + label metadata
│   ├── vespucci/                  # Vespucci results per dataset
│   ├── vespucci_leave_one_out/    # leave-one-replicate-out results (Calcagno2022)
│   ├── raw_data/                  # original raw data (e.g. Kathe2022 lamina annotations)
│   ├── rctd/                      # RCTD deconvolution results
│   │   ├── Calcagno2022/
│   │   └── regen_final/           # RCTD results for regeneration dataset
│   ├── regen/
│   │   ├── augur/                 # Augur feature importance results
│   │   └── aggregate_classifier/
│   ├── registration/              # registration robustness (real data)
│   │   ├── meta/
│   │   └── vespucci/
│   ├── spatial_cluster_genes/     # spatial gene-cluster analysis (regen)
│   │   └── gene_set_enrichment/
│   ├── summaries/
│   │   └── different_seeds/       # reproducibility-across-seeds summaries
│   └── DE_summaries/
│       ├── vespucci/
│       ├── vespucci_GO/
│       ├── others/
│       └── vespucci_de_genes_summary.rds  # cross-dataset DE gene summary (main Fig 5)
├── GO/
│   └── go.obo                     # Gene Ontology OBO file
├── metadata/
│   └── go_names.rds
├── kinematics/
│   └── plot_data.rds
└── Koupourtidou2024/              # external TBI bulk-RNA reference files (main Fig 2)
    ├── Belgard_2011_cortical_layer_genes.xls
    ├── Mus_musculus.NCBIM37.57.gtf.gz
    └── TBI_against_bulk_aucc_df.rds
```

### 3. R packages

Install dependencies from CRAN / Bioconductor before running any script:

```r
install.packages(c(
  'tidyverse', 'magrittr', 'ggplot2', 'patchwork', 'ggrastr', 'ggh4x',
  'pROC', 'PRROC', 'cetcolor', 'pals', 'paletteer', 'lawstat', 'nparcomp',
  'ontologyIndex', 'Matrix', 'scales'
))
BiocManager::install('Seurat')
```

---

## Reproducing each Figure

All scripts are in `R/figures/`. Run them from `~/git/vespucci-analysis/` or from any working directory (the `setwd` at the top of each script handles the path). Output PDFs are written to `fig/` inside this repo.

### Main Figures

| Figure | Script | Description |
|--------|--------|-------------|
| Fig 1 | `main-fig1.R` | Simulation overview: AUC map, runtime vs. Magellan, AUPR boxplot, delta-AUPR heatmap |
| Fig 2 | `main-fig2.R` | Koupourtidou2024: AUC map, top genes, GO modules, bulk-concordance lollipop |
| Fig 3 | `main-fig3.R` | Calcagno2022: registration, AUC map, top genes and GO modules |
| Fig 4 | `main-fig4.R` | Regeneration dataset: housekeeper genes, AUC maps, immune-cell RCTD, DE genes and GO modules per comparison |
| Fig 5 | `main-fig5.R` | Regeneration young-vs-old: spatial maps of top DE genes and per-method gene ranks |

### Extended Figures

| Figure | Script(s) | Description |
|--------|-----------|-------------|
| EFig 1 | `supp-simulations-example-genes.R` | Example spatial gene-expression maps coloured by TP/TN/FP/FN classification for each DE method |
| EFig 2 | `supp-simulations-performance.R` | Vespucci vs. Magellan run-time and region-detection accuracy; DE gene AUPRC and delta-AUPRC heatmap per simulation |
| EFig 3 | `supp-simulations-additional-metrics.R` | Accuracy, sensitivity and specificity of DE gene detection across simulations and methods |
| EFig 4 | `supp-simulations-false-discoveries.R` | False-discovery rate under a null simulation with no spatial signal |
| EFig 5 | `supp-simulations-resolutions-size.R` | Effect of spatial resolution and perturbed-region size on Vespucci performance |
| EFig 6 | `supp-simulations-tests.R` + `supp-real-data-tests.R` | Reproducibility across random seeds, spatial resolutions, subsampling depths, and registration shifts — both simulations and real datasets |
| EFig 7 | `supp-koupourtidou.R` | Koupourtidou2024: spatial expression of top DE genes and GO modules; lollipop plots |
| EFig 8 | `supp-calcagno.R` | Calcagno2022: RCTD cell-type annotation, spatial gene/GO-module maps, lollipop plots, leave-one-replicate-out robustness |
| EFig 9 | `supp-maniatis.R` | Maniatis2019: spatial AUC, top DE genes and GO modules, lollipop plots |
| EFig 10 | `supp-kathe.R` | Kathe2022: lamina annotation, spatial AUC, top DE genes and GO modules, lollipop plots |
| EFig 11 | `supp-zeng.R` | Zeng2023: registration, spatial AUC, gene maps, cell-type AUC comparison |
| EFig 12 | — | — |
| EFig 13 | `supp-regen.R` | Regeneration dataset: kinematics, spatial AUC maps, DE genes and GO modules across comparisons |
| EFig 14 | `supp-spatial-cluster-genes.R` | Spatially clustered gene-module analysis for the regeneration dataset |
| EFig 15 | `supp-registration-effect.R` | Registration-shift robustness experiments on simulated and real data |
| EFig 16 | `supp-feature-importance.R` | Feature-importance comparison between Vespucci, Augur and an aggregate classifier for the regeneration dataset |

---

## Notes

- Scripts load data via `data/`, which must be a symlink to the extracted data archive. Do not rename subdirectories.
- Output figures are saved to `fig/<Fig#>/` inside this repo. Create these directories before running if they do not exist.
- Scripts call `source('R/theme.R')` and (most) `source('R/functions/utils.R')` from this repo's `R/` directory.
