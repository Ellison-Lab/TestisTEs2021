library(tidyverse)
library(arrow)
library(ragg)
library(ComplexHeatmap)

# source("workflow/fig-scripts/theme.R")

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv')
w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')
w1118.scaled <- open_dataset("results/finalized/larval-w1118-testes/scaled", format='arrow')
te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

obstmp <- collect(w1118.obs) %>%
  left_join(rename.table) %>%
  mutate(clusters2 = clusters.rename)

clusters.nums <- obstmp %>% pull(clusters) %>% unique() %>% as.numeric %>% as.list %>% set_names(.,.)

te_expr_df <- map_df(clusters.nums,
                     ~{filter(w1118.scaled, (clusters == .) & (gene_id %in% unique(te.lookup$merged_te))) %>% collect()}) %>%
  dplyr::select(index, gene_id, expression, clusters) %>%
  mutate(clusters = as.character(clusters))

te_expr_mat <- te_expr_df %>% 
  dplyr::select(-clusters) %>%
  group_by(gene_id) %>%
  pivot_wider(names_from = index, values_from = expression) %>%
  column_to_rownames('gene_id') %>% 
  as.matrix()

palette_len <- 255
colors <- colorRampPalette( (c("#3B9AB2", "#3B9AB2","#78B7C5", "#E1AF00", "#F21A00")))(palette_len)
#colors <- gte21_pal("diverging")(palette_len)[c(1:10,40:50)]
#palette_len <- 21

myBreaks <- c(seq(min(te_expr_mat), mean(te_expr_mat), length.out=ceiling(palette_len/2)),
              seq(mean(te_expr_mat) + 0.01, 0.3*max(te_expr_mat), length.out=floor(palette_len/2)))

cls <- circlize::colorRamp2(breaks = myBreaks,colors = colors)

levs <- obstmp %>% 
  pull(clusters2) %>% unique %>%
  set_names(.,.) %>%
  map(~str_remove(.,'\\d+\\/')) %>%
  {names(.)[order(unlist(.))]} 

col <- list(cell_types=RColorBrewer::brewer.pal(length(levs),'Set3') %>%
              set_names(levs))

# Create the heatmap annotation
ha <- HeatmapAnnotation(cell_types = pull(obstmp,'clusters2'), col = col)


agg_png(snakemake@output[['png']], width=22, height =15, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
Heatmap(te_expr_mat[,obstmp$X1],
        column_names_rot = 45,
        column_title_side = "bottom",
        na_col = 'green',
        cluster_rows = T,
        name = 'Scaled expression',
        row_title = paste0('TE families (', nrow(te_expr_mat),")"), 
        #top_annotation = ha,
        show_column_names = F,
        col = cls,
        column_title_rot = 90,
        column_title_gp = gpar(fontsize=12),
        row_names_gp = gpar(fontsize = 2),show_row_names = F,
        show_column_dend = F, show_row_dend = F,
        column_split = factor(pull(obstmp,'clusters2'),levels=levs) %>% fct_reorder(.,as.numeric(str_extract(.,"\\d+"))),
        use_raster = T,
        heatmap_height = unit(4.7,units = "in"),
        width=unit(8.7,units="in"),
        cluster_column_slices = F)
dev.off()


te_expr_mat[,obstmp$X1] %>%
  as_tibble(rownames='TEs') %>%
  gather(cell,scaled.expression,-TEs) %>%
  write_tsv(snakemake@output[['dat']])

