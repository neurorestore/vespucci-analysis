setwd('~/git/vespucci-analysis/')
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(FactoMineR)
library(readxl)
source("R/theme.R")
source("R/kinematics/kinematics_parsing.R")

data_path = 'data/kinematics/final_data.rds'
dat = readRDS(data_path)
meta_cols = colnames(dat)[1:4]

p1_col = 'drag_duration_percentage'
p2_col = 'step_height_foot'
p3_col = 'ampl_elev_tight'

p1_name = "% Drag"
p2_name = "Step height"
p3_name = "Thigh oscillation"

p1_lab = "% Drag (% total time)"
p2_lab = "Step height of foot (cm)"
p3_lab = "Thigh oscillation (degrees)"

factor_colors = c(
    'young' = 'grey70',
    'old' = 'grey10',
    'treated' = nr_tree_red[70]
)

kin = parse_kinematics(dat,
                      # Meta data info
                      meta_cols = meta_cols,
                      label_col = 'condition',
                      replicate_col = 'replicate',
                      # Variables to plot
                      p1_col = p1_col,
                      p2_col = p2_col,
                      p3_col = p3_col,
                      pc_flip = F,
                      plot_data = T,
                      factor_colors = factor_colors)

p1 = plot_kinematics(kin, 
                     # Names of these variables for the loadings plot
                     p1_name = p1_name,
                     p2_name = p2_name,
                     p3_name = p3_name,
                     # Axis labels
                     p1_lab = p1_lab,
                     p2_lab = p2_lab,
                     p3_lab = p3_lab,
                     do_statistics = F,
                     factor_colors = factor_colors
)

p1 = p1 %>% wrap_plots(nrow = 1)
p1
