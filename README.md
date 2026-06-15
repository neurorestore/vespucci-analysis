# README

Scripts used to reproduce all figures in _"Identification of perturbation-responsive regions and genes in comparative spatial transcriptomics atlases"_.

## Repository layout

```
R/figures/          # one script per figure panel (see table below)
R/functions/        # shared utility functions
R/theme.R           # shared ggplot theme
```

## Reproducing figures

All scripts are in `R/figures/`. Each calls `setwd('~/git/vespucci-analysis/')` at the top and reads data via the `data/` symlink. Output PDFs are written to `fig/` inside this repo.

### Main Figures

| Figure | Script |
|--------|--------|
| Fig 1 | `main-fig1.R` |
| Fig 2 | `main-fig2.R` |
| Fig 3 | `main-fig3.R` |
| Fig 4 | `main-fig4.R` |
| Fig 5 | `main-fig5.R` |

### Extended Figures

| Figure | Script(s) |
|--------|-----------|
| EFig 1 | `supp-simulations-example-genes.R` |
| EFig 2 | `supp-simulations-performance.R` |
| EFig 3 | `supp-simulations-additional-metrics.R` |
| EFig 4 | `supp-simulations-false-discoveries.R` |
| EFig 5 | `supp-simulations-resolutions-size.R` |
| EFig 6 | `supp-simulations-tests.R` + `supp-real-data-tests.R` |
| EFig 7 | `supp-koupourtidou.R` |
| EFig 8 | `supp-calcagno.R` |
| EFig 9 | `supp-maniatis.R` |
| EFig 10 | `supp-kathe.R` |
| EFig 11 | `supp-zeng.R` |
| EFig 12 | — |
| EFig 13 | `supp-regen.R` |
| EFig 14 | `supp-spatial-cluster-genes.R` |
| EFig 15 | `supp-registration-effect.R` |
| EFig 16 | `supp-feature-importance.R` |
