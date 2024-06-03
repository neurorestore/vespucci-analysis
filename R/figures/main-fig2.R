setwd('~/git/vespucci-analysis/')
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(pROC)
library(PRROC)
library(lawstat)
library(nparcomp)
library(gridExtra)
library(cetcolor)
source('R/theme.R')

# housekeeper genes
sc = readRDS('data/regen/seurat/regen_final.rds')
sc@meta.data %<>%
    mutate(
        replicate = paste0(label, '_', slide)
    )

# extract metadata
meta = sc@meta.data %>%
    mutate(
        label = fct_recode(
            label,
            'Young' = 'young',
            'Old' = 'old',
            'Treated' = 'treated'
        ) %>%
            fct_relevel('Young', 'Old', 'Treated')
    )
expr = GetAssayData(sc, slot='counts') %>% NormalizeData()
housekeeper_genes = c('Gfap', 'Nefl', 'Mbp')
housekeeper_genes[!housekeeper_genes %in% rownames(expr)]

for (gene in housekeeper_genes) {
    if (gene == first(housekeeper_genes)) plot_list = list()
    meta$expr = expr[gene, meta$barcode]
    ## by condition
    conditions = unique(meta$label)
    coords = map_dfr(conditions, ~ {
        condition = .x
        fit = loess(expr ~ x * y, data = filter(meta, label == condition),
                    span = 0.02, degree = 1)
        coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5)) %>%
            mutate(interp = predict(fit, newdata = .)) %>%
            drop_na(interp) %>%
            mutate(condition = condition)
    })
    # range = range(coords$interp)
    range = quantile(coords$interp, c(0.01, 0.99))
    brks = c(range[1] + 0.1 * diff(range),
             range[2] - 0.1 * diff(range))
    labels = c('min', 'max')
    p1a = coords %>%
        mutate(
            interp = winsorize(interp, quantile(coords$interp, c(0.01, 0.99))),
            gene = gene
        ) %>%
        ggplot(aes(x = x, y = y, fill = interp)) +
        facet_grid(gene~ condition, switch='y') +
        geom_raster() +
        scale_y_reverse(expand = c(0, 0), breaks = c(-150, 150),
                        labels = c(expression('R'%->%""),
                                   expression(""%<-%'L')),
                        name = gene
        ) +
        scale_x_continuous(expand = c(0, 0), breaks = c(-400, 400),
                           labels = c(expression(""%<-%'Rostral'),
                                      expression('Caudal'%->%""))
        ) +
        scale_fill_gradientn(name = 'Expression', colors = nr_tree_red,
                             limits = range, breaks = brks, labels = labels) +
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
        coord_fixed() +
        boxed_theme() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(angle = 90, hjust = c(1, 0)),
            axis.text.x = element_text(hjust = c(0, 1)),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.length.x = unit(0, 'lines'),
            axis.ticks.length.y = unit(0, 'lines'),
            legend.position = 'right',
            legend.justification = 'bottom',
            legend.key.width = unit(0.18, 'lines'),
            legend.key.height = unit(0.18, 'lines')
        )
    if (gene != first(housekeeper_genes)) {
        p1a = p1a + theme(strip.text.x = element_blank())
    }
    if (gene != last(housekeeper_genes)) {
        p1a = p1a + theme(legend.position = 'none')
    }
    p1a
    plot_list[[length(plot_list)+1]] = p1a
}
p1 = wrap_plots(plot_list, ncol=1)
p1
ggsave('fig/Fig2/marker-genes.pdf', p1, width = 10.5, height = 5, units='cm')

ves_dir = 'data/regen/vespucci/'
comparison_names = c(
    'young_old' = 'Young vs. old',
    'treated_old' = 'Treated vs. old'
)

fixed_coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5))
plot_list = list()
for (comparison_name in names(comparison_names)) {
    ves_file = paste0(ves_dir, 'regen_final_', comparison_name, '-seed=42-nsub=10.rds')
    ves_res = readRDS(ves_file)$spatial_auc_result$aucs

    dat0 = ves_res %>%
        left_join(meta %>%
                      mutate(barcode = gsub('-', '_', barcode)) %>%
                      dplyr::select(x, y, label, barcode)
        )
    fit = loess(auc ~ x * y, data = dat0, span = 0.02, degree = 1)

    dat1 = fixed_coords %>%
        mutate(
            auc_fit = predict(fit, .),
            compar = unname(comparison_names[comparison_name])
        ) %>%
        filter(!is.na(auc_fit))

    range = range(dat1$auc_fit)
    brks = c(range[1] + 0.1 * diff(range),
             range[2] - 0.1 * diff(range))
    labels = format(range, digits = 2)
    new_color_pal = cet_pal(100, name = 'l19') %>% rev()

    p2_1 = dat1 %>%
        arrange(auc_fit) %>%
        ggplot(aes(x = x, y = y, fill = auc_fit, color = auc_fit)) +
        ggrastr::rasterise(
            geom_point(size = 0.5, shape = 21, stroke = 0), dpi = 300
        ) +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        scale_fill_gradientn(colours = new_color_pal,
                             name = 'AUC', labels = labels,
                             limits = range, breaks = brks) +
        scale_color_gradientn(colours = new_color_pal,
                              name = 'AUC', labels = labels,
                              limits = range, breaks = brks) +
        guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE),
               color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
        coord_fixed() +
        boxed_theme() +
        facet_wrap(~compar, nrow=1) +
        # geom_vline(xintercept = 0, size=0.1, linetype='dashed') +
        # geom_hline(yintercept = 0, size=0.1, linetype='dashed') +
        theme(
            # aspect.ratio = 0.5,
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            # axis.text.y = element_text(angle = 90, hjust = c(1, 0)),
            # axis.text.x = element_text(hjust = c(0, 1)),
            axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.length.x = unit(0, 'lines'),
            axis.ticks.length.y = unit(0, 'lines'),
            legend.position = 'right',
            legend.justification = 'bottom',
            legend.key.width = unit(0.18, 'lines'),
            legend.key.height = unit(0.18, 'lines'),
            plot.title = element_text(size = 5)
        )
    plot_list[[length(plot_list)+1]] = p2_1
}
p2 = wrap_plots(plot_list, ncol=2)
# p2
ggsave('fig/Fig2/auc.pdf', p2, width = 8, height = 6, units='cm')

rctd_res = readRDS('data/regen/analysis/summaries/regen-rctd-res-layer1.rds') %>%
    filter(
        grepl('treated', sp) | grepl('v14', sp)
    ) %>%
    mutate(
        label = ifelse(ref != '7d', 'Old', 'Young'),
        label = ifelse(grepl('treated', sp), 'Treated', label),
        barcode = gsub('-', '_', barcode)
    ) %>%
    left_join(meta %>%
                  mutate(barcode = gsub('-', '_', barcode)) %>%
                  dplyr::select(barcode, x, y, label)
              ) %>%
    mutate(label = factor(label, levels=c('Young', 'Old', 'Treated'))) %>%
    filter(
        !is.na(x),
        !is.na(y)
    ) %>%
    as.data.frame()

curr_color_pal['Immune cells'] = 'black'

rctd_res$cell_type = factor(rctd_res$cell_type, levels=names(curr_color_pal))
rctd_res$cell_type_bin = ifelse(as.character(rctd_res$cell_type == 'Immune cells'), 1, 0)

fixed_coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5))
grad = colorRampPalette(c('grey88', 'grey98', curr_color_pal[['Immune cells']]))(3)

plot_df = data.frame()
for (label in levels(rctd_res$label)) {
    rctd_res1 = rctd_res %>%
        filter(
            label == !!label
        )
    # rctd_res1$weight = -log(rctd_res1$weight)
    fit = loess(cell_type_bin ~ x * y, data = rctd_res1, span = 0.02, degree = 1)
    plot_df %<>% rbind(
        fixed_coords %>%
            # rctd_res1 %>%
            mutate(
                cell_type_bin_fitted = predict(fit, .),
                # cell_type_bin_fitted = scale(cell_type_bin_fitted)[,1],
                label = label
            )
    )
}

plot_df$cell_type_bin_fitted = as.numeric(plot_df$cell_type_bin_fitted)
plot_df %<>% filter(!is.na(cell_type_bin_fitted))

range = range(plot_df$cell_type_bin_fitted)
brks = c(range[1] + 0.1 * diff(range), range[2] - 0.1 * diff(range))
labels = c('min', 'max')
plot_df$label = factor(plot_df$label, levels = levels(rctd_res$label))
p3 = plot_df %>%
    ggplot(aes(x = x, y = y, fill = cell_type_bin_fitted, color = cell_type_bin_fitted)) +
    ggrastr::rasterise(
        geom_point(size = 0.4, shape = 21, stroke = 0), dpi = 300
    ) +
    scale_y_reverse(expand = c(0, 0), breaks = c(-150, 150)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(-400, 400)) +
    scale_fill_gradientn(colours = grad,
                         name = 'Weight',
                         labels = labels,
                         limits = range,
                         breaks = brks
    ) +
    scale_color_gradientn(colours = grad,
                          name = 'Weight',
                          labels = labels,
                          limits = range,
                          breaks = brks
    ) +
    guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE),
           color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
    coord_fixed() +
    boxed_theme() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = c(1, 0)),
        axis.text.x = element_text(hjust = c(0, 1)),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0, 'lines'),
        axis.ticks.length.y = unit(0, 'lines'),
        legend.position = 'bottom',
        legend.justification = 'right',
        legend.key.width = unit(0.18, 'lines'),
        legend.key.height = unit(0.18, 'lines'),
        plot.title = element_text(size = 5)) +
    facet_wrap(~label, nrow=1)
p3
ggsave(paste0('fig/Fig2/immune-cells-rctd.pdf'), p3, width=8, height=5, units='cm')

expr = GetAssayData(sc, slot='counts') %>% NormalizeData()
meta = sc@meta.data %>%
    mutate(
        label = fct_recode(
            label,
            'Young' = 'young',
            'Old' = 'old',
            'Treated' = 'treated'
        ) %>% as.character()
    )

go_df = readRDS('data/metadata/go_names.rds') %>%
    mutate(go_sub = gsub('\\:', '-', go))

comparisons = list(
    'young-old'=list(
        go_suffix = 'young_old',
        factor_lvls = c('Young', 'Old'),
        genes = c(
            'C1qa',
            'Ndrg2',
            'Ctsd',
            'Ftl1'
        ),
        gos = c(
            'inflammatory response',
            'axonogenesis'
        )
    ),
    'treated-old'=list(
        go_suffix = 'treated_old',
        factor_lvls = c('Old', 'Treated'),
        genes = c(
            'Spp1',
            'Ndrg2',
            'Ctsd',
            'Ftl1'
        ),
        gos = c(
            'inflammatory response',
            'axonogenesis'
        )
    )
)

seurat_go_dir = 'data/regen/seurat_GO/'
for (comparison in names(comparisons)) {
    genes = comparisons[[comparison]]$genes
    meta0 = meta %>%
        filter(label %in% comparisons[[comparison]]$factor_lvls) %>%
        mutate(label = fct_relevel(label, comparisons[[comparison]]$factor_lvls))
    for (gene in genes) {
        if (gene == first(genes)) plot_list = list()
        meta0$expr = expr[gene, meta0$barcode]
        ## by condition
        conditions = levels(meta0$label)
        coords = map_dfr(conditions, ~ {
            condition = .x
            fit = loess(expr ~ x * y, data = filter(meta0, label == condition),
                        span = 0.02, degree = 1)
            coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5)) %>%
                mutate(interp = predict(fit, newdata = .)) %>%
                drop_na(interp) %>%
                mutate(condition = condition)
        })
        # range = range(coords$interp)
        range = quantile(coords$interp, c(0.01, 0.99))
        brks = c(range[1] + 0.1 * diff(range),
                 range[2] - 0.1 * diff(range))
        labels = c('min', 'max')
        pal = nr_heat_red_no_white %>% tail(-5)
        p4a = coords %>%
            mutate(
                interp = winsorize(interp, quantile(coords$interp, c(0.01, 0.99))),
                gene = gene
            ) %>%
            ggplot(aes(x = x, y = y, fill = interp)) +
            facet_grid(~condition) +
            geom_raster() +
            scale_y_reverse(expand = c(0, 0), breaks = c(-150, 150),
                            labels = c(expression('R'%->%""),
                                       expression(""%<-%'L')),
                            name = gene
            ) +
            scale_x_continuous(expand = c(0, 0), breaks = c(-400, 400),
                               labels = c(expression(""%<-%'Rostral'),
                                          expression('Caudal'%->%""))
            ) +
            scale_fill_gradientn(name = 'Expression', colors = pal,
                                 limits = range, breaks = brks, labels = labels) +
            guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
            coord_fixed() +
            ggtitle(gene) +
            boxed_theme() +
            theme(
                plot.title = element_text(size = 6, margin = margin(0,0,-5,0)),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_text(angle = 90, hjust = c(1, 0)),
                axis.text.x = element_text(hjust = c(0, 1)),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.length.x = unit(0, 'lines'),
                axis.ticks.length.y = unit(0, 'lines'),
                legend.position = 'none',
                legend.justification = 'bottom',
                legend.key.width = unit(0.18, 'lines'),
                legend.key.height = unit(0.18, 'lines')
            )
        p4a
        plot_list[[length(plot_list)+1]] = p4a
    }

    go_sc = readRDS(paste0(seurat_go_dir, 'regen_final_', comparisons[[comparison]]$go_suffix, '.rds'))
    go_expr = GetAssayData(go_sc)
    gos = comparisons[[comparison]]$gos

    meta0 %<>% mutate(barcode = gsub('-', '_', barcode))

    for (go in gos) {
        go_sub = go_df %>%
            filter(go_name == !!go) %>%
            pull(go_sub)

        meta0$expr = go_expr[go_sub, meta0$barcode]
        ## by condition
        conditions = levels(meta0$label)
        coords = map_dfr(conditions, ~ {
            condition = .x
            fit = loess(expr ~ x * y, data = filter(meta0, label == condition),
                        span = 0.02, degree = 1)
            coords = tidyr::crossing(x = seq(-400, 400, 5), y = seq(-150, 150, 5)) %>%
                mutate(interp = predict(fit, newdata = .)) %>%
                drop_na(interp) %>%
                mutate(condition = condition)
        })
        # range = range(coords$interp)
        range = quantile(coords$interp, c(0.01, 0.99))
        brks = c(range[1] + 0.1 * diff(range),
                 range[2] - 0.1 * diff(range))
        labels = c('min', 'max')
        pal = nr_heat_blue_spatial
        p4a = coords %>%
            mutate(
                interp = winsorize(interp, quantile(coords$interp, c(0.01, 0.99))),
                gene = gene
            ) %>%
            ggplot(aes(x = x, y = y, fill = interp)) +
            facet_grid(~condition) +
            geom_raster() +
            scale_y_reverse(expand = c(0, 0), breaks = c(-150, 150),
                            labels = c(expression('R'%->%""),
                                       expression(""%<-%'L')),
                            name = gene
            ) +
            scale_x_continuous(expand = c(0, 0), breaks = c(-400, 400),
                               labels = c(expression(""%<-%'Rostral'),
                                          expression('Caudal'%->%""))
            ) +
            scale_fill_gradientn(name = 'Expression', colors = pal,
                                 limits = range, breaks = brks, labels = labels) +
            guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
            coord_fixed() +
            ggtitle(go) +
            boxed_theme() +
            theme(
                plot.title = element_text(size = 6, margin = margin(0,0,-5,0)),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_text(angle = 90, hjust = c(1, 0)),
                axis.text.x = element_text(hjust = c(0, 1)),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.length.x = unit(0, 'lines'),
                axis.ticks.length.y = unit(0, 'lines'),
                legend.position = 'none',
                legend.justification = 'bottom',
                legend.key.width = unit(0.18, 'lines'),
                legend.key.height = unit(0.18, 'lines')
            )
        p4a
        plot_list[[length(plot_list)+1]] = p4a
    }

    p4 = wrap_plots(plot_list, nrow=2)
    p4
    ggsave(paste0('fig/Fig2/', comparison, '-genes-gos.pdf'), p4, width = 17.5, height = 4.9, units='cm')
}
