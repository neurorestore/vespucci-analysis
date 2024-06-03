# Functions to process kinematics data quickly

# Function to parse kinematics data, and return an object containing
## Outcome measures matrix
## PC scores
## Meta data
library(FactoMineR)
parse_kinematics = function(input, 
                            # Define the columns to be split off that contain
                            # meta data
                            meta_cols = NULL,
                            label_col = NULL,
                            replicate_col = NULL,
                            # Define the columns to include in the figure
                            # right now just three can be picked
                            p1_col = NULL,
                            p2_col = NULL,
                            p3_col = NULL,
                            # Do we want to flip the PC direction?
                            pc_flip = F,
                            # Do we want to show the data?
                            plot_data = F,
                            # Do we need to drop any labels out?
                            drop_labels = NULL, 
                            factor_colors = NULL) {
  # split off the meta data and matrix
  meta = input %>% dplyr::select(all_of(meta_cols)) %>%
    dplyr::rename(label = all_of(label_col),
                  replicate = all_of(replicate_col)) %>%
    mutate(replicate = paste0(label, '-', replicate))
  mat = input %>% dplyr::select(-all_of(meta_cols)) %>%
    dplyr::rename(p1_col = all_of(p1_col),
                  p2_col = all_of(p2_col),
                  p3_col = all_of(p3_col))
  
  # optionally, drop some labels
  if (!is.null(drop_labels)) {
    keep = !meta$label %in% drop_labels
    meta = meta[keep, ]
    mat = mat[keep,]
  }

  # run PCA
  pca = FactoMineR::PCA(mat, graph = F)
  
  # parse PCA data
  pdat = pca$ind$coord %>% as.data.frame() %>%
    dplyr::select(1:2) %>%
    set_colnames(c("PC1", "PC2")) %>%
    ungroup()
  
  # optionally, flip the PC1 direction
  if (pc_flip) {
    pdat %<>% mutate(PC1 = -PC1)
  }
  
  # optionally, plot the data
  if (plot_data) {
    p1 = ggplot(pdat %>% bind_cols(meta), 
           aes(x = PC1, y = PC2, color = label)) +
      geom_point() +
      boxed_theme()
    if (!is.null(factor_colors)) {
        p1 = p1 + scale_color_manual(values = factor_colors)
    }
    
    print(p1)
  }
  
  # define the correlation of each variable to PC1/PC2
  loadings = pca$var$coord %>% as.data.frame() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'variable') %>%
    dplyr::rename(PC1 = "Dim.1", PC2 = "Dim.2") %>%
    dplyr::select(variable, PC1, PC2) %>%
    # dummy variable for X
    mutate(x = 'PC1') %>%
    mutate(PC1 = abs(PC1)) %>%
    arrange(desc(PC1))
  
  # get the variance explained by PC1 and PC2 for plotting
  var_expl = pca$eig %>% as.data.frame() %>%
    rownames_to_column(var = 'PC') %>%
    mutate(PC = gsub("comp ", "PC", PC)) %>%
    dplyr::select(PC, `percentage of variance`) %>%
    dplyr::rename(var = `percentage of variance`)
  
  # output object
  obj = list(
    mat = mat,
    meta = meta,
    pcs = pdat,
    loadings = loadings,
    var = var_expl
  )
  return(obj)
}

# function to paste stat snippets
paste_stat_snippet = function(results) {
  x = paste(rownames(results), '=', signif(results[,4], digits = 3))
  x = paste(x, collapse = '; ')
  message(x)
}

# Function to create our standard set of 6 plots
# define the palette
plot_kinematics = function(
    input,
    cols = nr_base_4,
    # Define the names for the loading plot
    p1_name = NULL,
    p2_name = NULL,
    p3_name = NULL,
    # Define the axis labels
    p1_lab = NULL,
    p2_lab = NULL,
    p3_lab = NULL,
    PCA_y_label = 'Walking score (Scaled PC1)',
    factor_colors = NULL,
    do_statistics = T
) {
  
  if (is.null(factor_colors)) {
      levels = unique(input$condition)
  } else {
      levels = names(factor_colors)
  }
  
  # summarize mean of PC1/PC2 for each group
  pcs0 = input$pcs %>%
    bind_cols(input$meta %>% dplyr::select(label, replicate)) %>%
    group_by(label, replicate) %>%
    summarise_at(vars(PC1, PC2), mean) %>%
    # rescale these
    ungroup() %>%
    mutate_at(vars(PC1, PC2), function(x) scales::rescale(x, to = c(0,1)))
  
  if (!is.null(levels)) {
    pcs0$label = factor(pcs0$label, levels=levels)
  }
  
  sum = pcs0 %>%
    group_by(label) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2))
  
  # prepare the outcomes for plotting
  outcomes = input$mat %>%
    bind_cols(input$meta %>% dplyr::select(label, replicate)) %>%
    group_by(label, replicate) %>%
    summarise_all(mean)
  
  # define the PC1/PC2 labels
  x_lab = paste0("Scaled PC1 (", input$var %>% 
                   dplyr::filter(PC == 'PC1') %>% 
                   pull(var) %>%
                   round(digits = 2),
                 '%)')
  y_lab = paste0("Scaled PC2 (", input$var %>% 
                   dplyr::filter(PC == 'PC2') %>% 
                   pull(var) %>%
                   round(digits = 2),
                 '%)')
  
  if (do_statistics) {
    if (is.null(levels)) levels = c("uninjured", "acute", "chronic")
    # get the factor levels
    label_f = factor(pcs0$label, levels = levels)
    fct_levels = crossing(
      start = levels(label_f),
      end = levels(label_f)
    ) %>%
      left_join(
        data.frame(
          start = levels(label_f),
          start_x = seq_len(length(levels(label_f)))
        )
      ) %>%
      left_join(
        data.frame(
          end = levels(label_f),
          end_x = seq_len(length(levels(label_f)))
        )
      )
  }
  
  # First plot, loading ranking
  ### loadings
  loadings = input$loadings %>%
    mutate(variable = str_replace(variable, 'p1_col', p1_name)) %>%
    mutate(variable = str_replace(variable, 'p2_col', p2_name)) %>%
    mutate(variable = str_replace(variable, 'p3_col', p3_name)) %>%
    mutate(idx = row_number(),
           color = variable %in% c(p1_name, p2_name, p3_name))
  labs = loadings %>% 
    dplyr::filter(variable %in% c(p1_name, p2_name, p3_name))
  
  p1 = loadings %>% 
    ggplot(aes(x = idx, y = PC1, color = color)) +
    geom_segment(aes(xend = idx, y = 0, yend = PC1), size = 0.2) +
    geom_text_repel(data = labs, aes(label = variable), size = 1.75,
                    # label.size = NA, label.padding = unit(0.4, 'lines'),
                    color = 'black', segment.size = 0.2,
                    min.segment.length = 0.1, nudge_y = .35) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_reverse(breaks = c(1, 25, 50, 75, 100)) +
    guides(color = 'none', fill = 'none') +
    labs(x = 'Variable rank', y = 'Correlation to PC1') +
    scale_fill_manual('', values = c('grey60', 'indianred3'), labels = parse_format()) +
    scale_color_manual('', values = c('grey60', 'indianred3'), labels = parse_format()) +
    coord_cartesian(ylim = c(0, 1.5)) +
    boxed_theme(size_lg = 5, size_sm = 5) +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          panel.border = element_rect(size = 0.2),
          axis.ticks.x = element_line(size = 0.2),
          axis.ticks.y = element_line(size = 0.2),
          legend.position = 'none',
          legend.key.size = unit(0.6, "lines"))
  p1
  
  # next plot, showing the PCs as a dot plot
  # pcs0$label = factor(pcs0$label, levels = c("uninjured", "acute", "chronic" ))
  
  pcs0$label = factor(pcs0$label, levels = levels)
  p2 = ggplot(pcs0, aes(x = PC1, y = PC2, color = label, fill = label)) +
    geom_point(size = 0.3, shape = 21, stroke = 0) +
    geom_point(data = sum, size = 2, show.legend = F, shape = 21, 
               stroke = 0.15, color = 'grey20') +
    # scale_fill_manual('', values = c('grey70', 'grey55', 'grey34', 'grey18' , 'indianred3'), labels = parse_format()) +
    # scale_color_manual('', values = c('grey70', 'grey55', 'grey34', 'grey18' , 'indianred3'), labels = parse_format()) +
    scale_fill_manual(values = factor_colors) +
    scale_color_manual(values = factor_colors) +
    scale_y_continuous(y_lab) +
    scale_x_continuous(x_lab) +
    # coord_fixed() +
    guides(color = guide_legend(override.aes = list(size = 0.7, alpha = 1))) +
    grid_theme(size_lg = 5, size_sm = 5) +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          panel.border = element_rect(size = 0.2),
          axis.ticks.x = element_line(size = 0.2),
          axis.ticks.y = element_line(size = 0.2),
          legend.position = 'none',
          legend.key.size = unit(0.6, "lines"))
  
  p2
  
  # plot the PCs as a bar graph now
  max_pc = round(max(pcs0$PC1))
  # min_pc = round(min(pcs0$PC1) - 2)
  min_pc = 0
  
  p3 = pcs0 %>% 
    dplyr::select(label, replicate, PC1) %>%
    group_by(label, replicate) %>%
    dplyr::rename(variable = PC1) %>%
    summarise(mean = mean(variable)) %>%
    ggplot(aes(x = label, y = mean, color = label, fill = label)) +
    stat_summary(geom = "bar", fun = mean, alpha = 0.7, width = 0.7, 
                 size = 0.3, color = NA) +
    geom_hline(aes(yintercept = 0), color = 'grey50', size = 0.2) +
    geom_jitter(size = 0.5, width = 0.05, height = 0) +
    # scale_fill_manual('', values = c('grey70', 'grey55', 'grey34', 'grey18' , 'indianred3'), labels = parse_format()) +
    # scale_color_manual('', values = c('grey70', 'grey55', 'grey34', 'grey18' , 'indianred3'), labels = parse_format()) +
    scale_fill_manual(values = factor_colors) +
    scale_color_manual(values = factor_colors) +
    scale_y_continuous(PCA_y_label,
                       expand = c(0, 0),
                       # limits = c(min_pc, max_pc),
                       # breaks = seq(min_pc, max_pc, 2)
                       limits = c(0, 1),
                       breaks = seq(0, 1, 1)
                       ) +
    guides(fill = 'none') +
    clean_theme(size_lg = 5, size_sm = 5) +
    theme(legend.position = 'none',
      aspect.ratio = 2.5,
      plot.background = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_line(size = 0.2, color = 'grey50'),
      axis.ticks.y = element_line(size = 0.2),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  p3
  
  if (do_statistics) {
    # get the statistics
    pc_stat = aov(lm(PC1 ~ label, data = pcs0)) %>% TukeyHSD()
    paste_stat_snippet(pc_stat$label)
    stat_df = pc_stat$label %>%
      as.data.frame() %>%
      rownames_to_column(var = 'comparison') %>%
      separate(comparison, c('start', 'end'), '-') %>%
      mutate(sig = ifelse(`p adj` < 0.001, '***', 
                          ifelse(`p adj` < 0.01, '**',
                                 ifelse(`p adj` < 0.05, '*',
                                        'n.s.')))) %>%
      left_join(fct_levels) %>%
      filter(sig != 'n.s.')
    
    # add each stat in turn
    y_start = max(pcs0$PC1) + 0.3 * max(pcs0$PC1)
    map(seq_len(nrow(stat_df)), ~ {
      idx = .
      start_x = stat_df[idx,]$start_x
      end_x = stat_df[idx,]$end_x
      sig = stat_df[idx,]$sig
      y = y_start + ((0.3 * max(pcs0$PC1)) * (idx-1))
      y_end = y + (0.3 * max(pcs0$PC1))
      nudge = (0.1 * max(pcs0$PC1))
      stat_positions = seq(1.5, 5, 0.5)
      stat_pos = stat_positions[idx]
      p3 <<- p3 + 
        geom_segment(data = pcs0[1,],
                     aes(x = start_x, y = y, xend = end_x, yend = y), 
                     color = 'black',
                     size = 0.2) + 
        geom_text(data = pcs0[1,], aes(x = stat_pos, y = y, label = sig),
                  color = 'black', size = 1.75, nudge_y = nudge) +
        scale_y_continuous(PCA_y_label,
                           expand = c(0, 0),
                           limits = c(min_pc, y_end),
                           breaks = seq(min_pc, y_end, 4))
      })
      p3
    }

  ### First outcome measure
  p1_dat = outcomes %>% 
    ungroup() %>%
    filter(!is.na(p1_col))
    #mutate(p1_col = ifelse(p1_col < 1, max(p1_col, na.rm = T) * 0.01, p1_col))
  # p1_dat$label = factor(p1_dat$label, levels = c("uninjured", "acute", "chronic" ))
  if (!is.null(levels)) p1_dat$label = factor(p1_dat$label, levels = levels)
  p4 = p1_dat %>%
    ggplot(aes(x = label, y = p1_col, color = label, fill = label)) +
    stat_summary(geom = "bar", fun = mean, alpha = 0.7, width = 0.7, 
                 color = NA) +
    geom_jitter(size = 0.5, width = 0.05, height = 0) +
    # scale_fill_manual('', values = c('grey70', 'grey55', 'grey34', 'grey18' , 'indianred3'), labels = parse_format()) +
    # scale_color_manual('', values = c('grey70', 'grey55', 'grey34', 'grey18' , 'indianred3'), labels = parse_format()) +
    scale_fill_manual(values = factor_colors) +
    scale_color_manual(values = factor_colors) +
    guides(fill = 'none') +
    labs(y = p1_lab) +
    clean_theme(size_lg = 5, size_sm = 5) +
    theme(legend.position = 'none',
      plot.background = element_blank(),
      axis.title.x = element_blank(),
      aspect.ratio = 2.5,
      axis.line.y = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.y = element_line(size = 0.2),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  p4
  
  if (do_statistics) {
    # get the statistics
    p1_stat = aov(lm(p1_col ~ label, data = p1_dat)) %>% TukeyHSD()
    paste_stat_snippet(p1_stat$label)
    p1_stat_df = p1_stat$label %>%
      as.data.frame() %>%
      rownames_to_column(var = 'comparison') %>%
      separate(comparison, c('start', 'end'), '-') %>%
      mutate(sig = ifelse(`p adj` < 0.001, '***', 
                          ifelse(`p adj` < 0.01, '**',
                                 ifelse(`p adj` < 0.05, '*',
                                        'n.s.')))) %>%
      left_join(fct_levels) %>%
      filter(sig != 'n.s.')
    
    # add each stat in turn
    y_start = max(p1_dat$p1_col) + 0.1 * max(p1_dat$p1_col)
    map(seq_len(nrow(p1_stat_df)), ~ {
      idx = .
      start_x = p1_stat_df[idx,]$start_x
      end_x = p1_stat_df[idx,]$end_x
      sig = p1_stat_df[idx,]$sig
      y = y_start + ((0.1 * max(p1_dat$p1_col)) * (idx-1))
      y_end = y + (0.1 * max(p1_dat$p1_col))
      nudge = (0.05 * max(p1_dat$p1_col))
      stat_positions = seq(1.5, 5, 0.5)
      stat_pos = stat_positions[idx]
      p4 <<- p4 + 
        geom_segment(data = p1_dat[1,],
                     aes(x = start_x, y = y, xend = end_x, yend = y), 
                     color = 'black',
                     size = 0.2) + 
        geom_text(data = p1_dat[1,], aes(x = stat_pos, y = y, label = sig),
                  color = 'black', size = 1.75, nudge_y = nudge)
    })
    p4
  }
  
  # fix the breaks
  p4_brks = layer_scales(p4)$y$break_positions()
  if (is.na(last(p4_brks))) {
    p4_last_brk = last(na.omit(p4_brks)) + diff(tail(na.omit(p4_brks, 2))) %>% unique()
    p4_brks %<>% na.omit() %>% c(p4_last_brk)
  }
  # make sure the last break is higher than the max
  if (max(p1_dat$p1_col) > last(p4_brks)) {
    p4_last_brk = last(na.omit(p4_brks)) + diff(tail(na.omit(p4_brks, 2))) %>% unique()
    p4_brks %<>% c(p4_last_brk)
  }
  p4 = p4 +
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = last(p4_brks)),
                 color = 'grey50', size = 0.2)
  p4
  
  ## Second outcomes plot
  p2_dat = outcomes %>% filter(!is.na(p2_col))
  if (!is.null(levels)) p2_dat$label = factor(p2_dat$label, levels = levels)
  # p2_dat$label = factor(p2_dat$label, levels = c("uninjured", "acute", "chronic" ))
  p5 = p2_dat %>%
    ggplot(aes(x = label, y = p2_col, color = label, fill = label)) +
    stat_summary(geom = "bar", fun = mean, alpha = 0.7, width = 0.7, 
                 size = 0.3, color = NA) +
    geom_jitter(size = 0.5, width = 0.05, height = 0) +
    # scale_fill_manual('', values = c('grey70', 'grey55', 'grey34', 'grey18' , 'indianred3'), labels = parse_format()) +
    # scale_color_manual('', values = c('grey70', 'grey55', 'grey34', 'grey18' , 'indianred3'), labels = parse_format()) +
    scale_fill_manual(values = factor_colors) +
    scale_color_manual(values = factor_colors) +
    guides(fill = 'none') +
    labs(y = p2_lab) +
    clean_theme(size_lg = 5, size_sm = 5) +
    theme(legend.position = 'none',
          plot.background = element_blank(),
          axis.title.x = element_blank(),
          aspect.ratio = 2.5,
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.y = element_line(size = 0.2),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  p5
  
  if (do_statistics) {
    # get the statistics
    p2_stat = aov(lm(p2_col ~ label, data = p2_dat)) %>% TukeyHSD()
    paste_stat_snippet(p2_stat$label)
    p2_stat_df = p2_stat$label %>%
      as.data.frame() %>%
      rownames_to_column(var = 'comparison') %>%
      separate(comparison, c('start', 'end'), '-') %>%
      mutate(sig = ifelse(`p adj` < 0.001, '***', 
                          ifelse(`p adj` < 0.01, '**',
                                 ifelse(`p adj` < 0.05, '*',
                                        'n.s.')))) %>%
      left_join(fct_levels) %>%
      filter(sig != 'n.s.')
    
    # add each stat in turn
    y_start = max(p2_dat$p2_col) + 0.1 * max(p2_dat$p2_col)
    map(seq_len(nrow(p2_stat_df)), ~ {
      idx = .
      start_x = p2_stat_df[idx,]$start_x
      end_x = p2_stat_df[idx,]$end_x
      sig = p2_stat_df[idx,]$sig
      y = y_start + ((0.1 * max(p2_dat$p2_col)) * (idx-1))
      y_end = y + (0.1 * max(p2_dat$p2_col))
      nudge = (0.05 * max(p2_dat$p2_col))
      stat_positions = seq(1.5, 5, 0.5)
      stat_pos = stat_positions[idx]
      p5 <<- p5 + 
        geom_segment(data = p2_dat[1,],
                     aes(x = start_x, y = y, xend = end_x, yend = y), 
                     color = 'black',
                     size = 0.2) + 
        geom_text(data = p2_dat[1,], aes(x = stat_pos, y = y, label = sig),
                  color = 'black', size = 1.75, nudge_y = nudge)
    })
    p5
  }
  
  # fix the breaks
  p5_brks = layer_scales(p5)$y$break_positions()
  if (is.na(last(p5_brks))) {
    p5_last_brk = last(na.omit(p5_brks)) + diff(tail(na.omit(p5_brks, 2))) %>% unique()
    p5_brks %<>% na.omit() %>% c(p5_last_brk)
  }
  # make sure the last break is higher than the max
  if (max(p2_dat$p2_col) > last(p5_brks)) {
    p5_last_brk = last(na.omit(p5_brks)) + diff(tail(na.omit(p5_brks, 2))) %>% unique()
    p5_brks %<>% c(p5_last_brk)
  }
  p5 = p5 +
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = last(p5_brks)),
                 color = 'grey50', size = 0.2)
  p5
  
  ### Plot the 3rd outcome measure
  p3_dat = outcomes %>% dplyr::filter(!is.na(p3_col))
  if (!is.null(levels)) p3_dat$label = factor(p3_dat$label, levels = levels)
  # p3_dat$label = factor(p3_dat$label, levels = c("uninjured", "acute", "chronic" ))
  p6 = p3_dat %>%
    ggplot(aes(x = label, y = p3_col, color = label, fill = label)) +
    stat_summary(geom = "bar", fun = mean, alpha = 0.7, width = 0.7, 
                 size = 0.3, color = NA) +
    geom_jitter(size = 0.5, width = 0.05, height = 0) +
    # scale_fill_manual('', values = c('grey70', 'grey55', 'grey34', 'grey18' , 'indianred3'), labels = parse_format()) +
    # scale_color_manual('', values = c('grey70', 'grey55', 'grey34', 'grey18' , 'indianred3'), labels = parse_format()) +
    scale_fill_manual(values = factor_colors) +
    scale_color_manual(values = factor_colors) +
    guides(fill = 'none') +
    labs(y = p3_lab) +
    clean_theme(size_lg = 5, size_sm = 5) +
    theme(legend.position = 'right',
          plot.background = element_blank(),
          axis.title.x = element_blank(),
          aspect.ratio = 2.5,
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.y = element_line(size = 0.2),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  p6
  
  if (do_statistics) {
    # get the statistics
    p3_stat = aov(lm(p3_col ~ label, data = p3_dat)) %>% TukeyHSD()
    paste_stat_snippet(p3_stat$label)
    p3_stat_df = p3_stat$label %>%
      as.data.frame() %>%
      rownames_to_column(var = 'comparison') %>%
      separate(comparison, c('start', 'end'), '-') %>%
      mutate(sig = ifelse(`p adj` < 0.001, '***', 
                          ifelse(`p adj` < 0.01, '**',
                                 ifelse(`p adj` < 0.05, '*',
                                        'n.s.')))) %>%
      left_join(fct_levels) %>%
      filter(sig != 'n.s.')
    
    # add each stat in turn
    y_start = max(p3_dat$p3_col) + 0.1 * max(p3_dat$p3_col)
    map(seq_len(nrow(p3_stat_df)), ~ {
      idx = .
      start_x = p3_stat_df[idx,]$start_x
      end_x = p3_stat_df[idx,]$end_x
      sig = p3_stat_df[idx,]$sig
      y = y_start + ((0.1 * max(p3_dat$p3_col)) * (idx-1))
      y_end = y + (0.1 * max(p3_dat$p3_col))
      nudge = (0.05 * max(p3_dat$p3_col))
      stat_positions = seq(1.5, 5, 0.5)
      stat_pos = stat_positions[idx]
      p6 <<- p6 + 
        geom_segment(data = p3_dat[1,],
                     aes(x = start_x, y = y, xend = end_x, yend = y), 
                     color = 'black',
                     size = 0.2) + 
        geom_text(data = p3_dat[1,], aes(x = stat_pos, y = y, label = sig),
                  color = 'black', size = 1.75, nudge_y = nudge)
    })
    p6
  }
  
  # fix the breaks
  p6_brks = layer_scales(p6)$y$break_positions()
  # fix the bug causing the last break to show up as NA
  if (is.na(last(p6_brks))) {
    p6_last_brk = last(na.omit(p6_brks)) + diff(tail(na.omit(p6_brks, 2))) %>% unique()
    p6_brks %<>% na.omit() %>% c(p6_last_brk)
  }
  # make sure the last break is higher than the max
  if (max(p3_dat$p3_col) > last(p6_brks)) {
    p6_last_brk = last(na.omit(p6_brks)) + diff(tail(na.omit(p6_brks, 2))) %>% unique()
    p6_brks %<>% c(p6_last_brk)
  }
  p6 = p6 +
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = last(p6_brks)),
                 color = 'grey50', size = 0.2)
  p6
  
  # return the plots as a list in case any need to be fixed
  plots = list(p1, p2, p3, p4, p5, p6)
  return(plots)
}

# A function to return the statistics info to put in the figure legends
return_statistics_exerpt = function(input,
                                    p1_name = "% Drag",
                                    p2_name = "Foot oscillation",
                                    p3_name = "Acceleration at swing onset") {
  # summarise PCs
  pcs0 = input$pcs %>%
    bind_cols(input$meta %>% dplyr::select(label, replicate)) %>%
    group_by(label, replicate) %>%
    summarise_at(vars(PC1, PC2), mean)
  
  # do statistic for PC
  pc_stat = pcs0 %>% group_by(label, replicate) %>% summarise(PC1 = mean(PC1))
  tuk = aov(lm(PC1 ~ label, data = pc_stat)) %>% TukeyHSD()
  pc_snippet = paste(rownames(tuk$label), '=', signif(tuk$label[,4], digits = 3))
  pc_snippet = paste(pc_snippet, collapse = '; ')

  # do statistics for our three main outcomes
  outcomes_stat = input$mat %>%
    bind_cols(input$meta) %>%
    gather(variable, value, -replicate, -label) %>%
    mutate(value = as.numeric(value)) %>%
    filter(variable %in% c("p1_col", "p2_col", "p3_col")) %>%
    group_by(label, replicate, variable) %>% 
    summarise(value = mean(value, na.rm = T)) %>%
    ungroup() %>%
    split(.$variable)
  
  outcome_aovs = map(outcomes_stat, function(x) 
    aov(lm(value ~ label, data = x)) %>% TukeyHSD())
  outcomes_snippet = map(outcome_aovs, ~ {
    tmp = .
    snp = paste(rownames(tmp$label), '=', signif(tmp$label[,4], digits = 3))
    snp = paste(snp, collapse = '; ')
    return(snp)
  })
  
  # return a simple vector
  output = paste0(
    "PC1: ", pc_snippet, '; ',
    p1_name, ": ", outcomes_snippet$p1_col, '; ',
    p2_name, ": ", outcomes_snippet$p2_col, '; ',
    p3_name, ": ", outcomes_snippet$p3_col
  )
  return(output)
}
  
