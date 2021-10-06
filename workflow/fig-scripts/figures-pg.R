library(tidyverse)
library(plotgardener)
library(png)
library(grid)
library(magick)
library(ggforce)
library(ggtext)

set.seed(1)
source("workflow/fig-scripts/theme.R")

f1_fl <- "~/Desktop/figure1.pdf"
f2_fl <- "~/Desktop/figure2.pdf"
f3_fl <- "~/Desktop/figure3.pdf"
f4_fl <- "~/Desktop/figure4.pdf"
f5_fl <- "~/Desktop/figure5.pdf"

f1_fl <- snakemake@output[["f1"]]
f2_fl <- snakemake@output[["f2"]]
f3_fl <- snakemake@output[["f3"]]
f4_fl <- snakemake@output[["f4"]]
f5_fl <- snakemake@output[["f5"]]


# ------------- figure1 -----------------------------------------------------

umap <- read_rds('results/figs/intro_larval_umap/intro_larval_umap.ggp.rds') +
  theme(axis.text = element_text(size=unit(7,"pt")), axis.title = element_text(size=unit(7,"pt"))) +
  guides(color=F)

markers <- read_rds('results/figs/larval_marker_expression/larval_marker_expression.ggp.rds') +
  theme(aspect.ratio = NULL, text = element_text(size=unit(6,"pt"))) + 
  theme(axis.text.y = element_text(face="italic")) +
  scale_fill_distiller(palette = 7, direction = 1)

n_per_clust <- read_tsv('results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.dat.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename) %>%
  group_by(clusters.rename) %>%
  tally() %>%
  ungroup() %>%
  ggplot(aes(clusters.rename,n, fill=clusters.rename)) +
  geom_col() +
  theme_gte21() +
  scale_fill_gte21() +
  theme(axis.text.x = element_text(angle=90, hjust=1)) + xlab("") +
  ylab("cells per cluster") +
  guides(fill=F) + 
  theme(aspect.ratio = NULL) + 
  scale_y_continuous(breaks = seq(0,3000,by=500))

cairo_pdf(f1_fl,width = 7.1, height = 3.95)

pageCreate(width = 180, height = 100, default.units = "mm")

plotText(label = "a", fontsize = 7, x = 5, y = 5, just = "left", default.units = "mm")
plotText(label = "b", fontsize = 7, x = 95, y = 5, just = "left", default.units = "mm")
plotText(label = "c", fontsize = 7, x = 95, y = 60, just = "left", default.units = "mm")

plotGG(umap,x = 0,y = 5,width = 95, height = 90, default.units = "mm")
plotGG(markers,x = 90,y = 5,width = 85, height = 60, default.units = "mm")
plotGG(n_per_clust,x = 95,y = 55,width = 80, height = 50, default.units = "mm")

pageGuideHide()

dev.off()

# ------------- figure2 -----------------------------------------------------


heat1_img <- image_read('results/figs/larval_te_heatmap/larval_te_heatmap.png', density = 300,) %>% 
  image_trim(fuzz = 3) %>%
  image_background(color = "none")

heat1 <- image_ggplot(heat1_img)

te_expression_by_clust <- read_rds("results/figs/te_expression_by_cluster/te_expression_by_cluster.ggp.rds") + 
  theme(axis.title.x = element_blank())

cairo_pdf(f2_fl,width = 7.1, height = 2.55,fallback_resolution = 300)

pageCreate(width = 180, height = 65, default.units = "mm")

plotText(label = "a", fontsize = 7, x = 5, y = 5, just = "left", default.units = "mm")
plotText(label = "b", fontsize = 7, x = 120, y = 5, just = "left", default.units = "mm")

plotGG(heat1,x = 10,y = -10,width = 100, height = 90, default.units = "mm")
plotGG(te_expression_by_clust,x = 120,y = 5,width = 60, height = 60, default.units = "mm")

pageGuideHide()

dev.off()


# ------------- figure3 -----------------------------------------------------


umap_2 <- read_rds('results/figs/larval_tep_usage_umap/larval_tep_usage_umap.ggp.rds') +
  theme(aspect.ratio = NULL) +
  scale_color_distiller(type="div",palette = 4, name=str_wrap("Module score",width = 1))
  #theme(legend.position = c(0.9,0.13)) + coord_fixed() +

barch <- read_rds('results/figs/larval_all_gep_te_barchart/larval_all_gep_te_barchart.ggp.rds')  + theme(aspect.ratio = NULL) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  #theme(legend.position = c(0.1,0.75), legend.title = element_blank(), legend.margin = margin(1,1,1,1), legend.text = element_text(size=rel(0.5)), legend.background = element_rect(color = NA)) +
  xlab("Module")

piech <- read_rds('results/figs/larval_tep_pie/larval_tep_pie.ggp.rds') +
  ggtitle("Module 27") +
  theme(text=element_text(size=unit(5,"pt"))) +
  theme(legend.position = "bottom", legend.key.size = unit(1,"line"), legend.margin = margin(t=1,unit="line")) +
  theme(plot.margin = margin(), plot.title = element_text(hjust=0.5, margin = margin(5,5,5,5))) +
  theme(legend.direction = "horizontal", legend.position = c(0.52,0)) +
  ggtitle("Module 27") +
  theme(text=element_text(size=unit(5,"pt")), legend.text = element_text(size=unit(5,"pt")))

cairo_pdf(f3_fl,width = 7.1, height = 3.95)

pageCreate(width = 180, height = 100, default.units = "mm")

plotText(label = "a", fontsize = 7, x = 5, y = 5, just = "left", default.units = "mm")
plotText(label = "c", fontsize = 7, x = 85, y = 5, just = "left", default.units = "mm")
plotText(label = "b", fontsize = 7, x = 5, y = 60, just = "left", default.units = "mm")

plotGG(barch,x = 10,y = 5,width =75, height = 50, default.units = "mm")
plotGG(piech,x = 15,y = 56,width = 60, height = 30, default.units = "mm")
plotGG(umap_2,x = 90,y = 5,width = 90, height = 90, default.units = "mm")

pageGuideHide()

dev.off()


# ------------- figure4 -----------------------------------------------------

y_genes <- read_rds('results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.ggp.rds') +
  theme(aspect.ratio = NULL, strip.text = element_blank()) +
  theme(legend.text = element_text(size=unit(7,"pt")))

fish_te_umis <- read_rds("results/figs/larval_fish_candidate_umis/larval_fish_candidate_umis.ggp.rds") +
  theme(strip.text = element_text(face="italic"))

#ph <- image_read("/media/mlawlor/T7/microscopy_figs/210215_accord2_calfluor610_eachm_quasar670_3p4-2.slices_1_1.representative.png") %>% image_ggplot()
ph <- image_read("resources/210215_accord2_calfluor610_eachm_quasar670_3p4-2.slices_1_1.representative.png") %>% image_ggplot()

ph <- ph + geom_ellipse(aes(x0=1650, y0=1950, a=1550,b=1950, angle=pi*0.9), linetype=2, color="white")

cairo_pdf(f4_fl,width = 7.1, height = 7.1)

pageCreate(width = 180, height = 180, default.units = "mm")

plotText(label = "a", fontsize = 7, x = 5, y = 5, just = "left", default.units = "mm")
plotText(label = "b", fontsize = 7, x = 95, y = 5, just = "left", default.units = "mm")
plotText(label = "c", fontsize = 7, x = 5, y = 95, just = "left", default.units = "mm")

plotGG(fish_te_umis,x = 0,y = 5,width = 95, height = 90, default.units = "mm")
plotGG(ph,x = 90,y = 5,width = 85, height = 60, default.units = "mm")
plotGG(y_genes,x = 5,y = 95,width = 170, height = 90, default.units = "mm")


pageGuideHide()

dev.off()

# ------------- figure5 -----------------------------------------------------

larrac <- read_rds('results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.ggp.rds') + 
  theme(aspect.ratio = NULL) +
  theme(axis.title.x=element_blank())

tidal <- read_rds('results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.ggp.rds') + 
  theme(aspect.ratio = NULL) + guides(fill=F)

w1118_copies_box <- read_rds('results/figs/w1118_y_linked_copies/w1118_y_linked_copies.1.ggp.rds') + 
  theme(aspect.ratio = NULL)

male_expr <- read_rds('results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.ggp.rds') +
  theme(aspect.ratio = NULL) + guides(fill=F) 

pirna <- read_rds("results/figs/larval_pirna_expression/larval_pirna_expression.ggp.rds") +  
  theme(aspect.ratio = NULL) +
  theme(axis.text.y = element_markdown())

cairo_pdf(f5_fl,width = 7.1, height = 4.33)

pageCreate(width = 180, height = 110, default.units = "mm")

plotText(label = "a", fontsize = 7, x = 5, y = 5, just = "left", default.units = "mm")
plotText(label = "b", fontsize = 7, x = 55, y = 5, just = "left", default.units = "mm")
plotText(label = "c", fontsize = 7, x = 5, y = 50, just = "left", default.units = "mm")
plotText(label = "d", fontsize = 7, x = 55, y = 50, just = "left", default.units = "mm")
plotText(label = "e", fontsize = 7, x = 110, y = 5, just = "left", default.units = "mm")

plotGG(larrac, x = 3, y = 7, width = 55, height=45, default.units = "mm")
plotGG(w1118_copies_box,x = 55, y = 7,  width = 55, height=45,default.units = "mm")
plotGG(male_expr,x = 3, y = 52, width = 55, height=45,default.units = "mm")
plotGG(tidal,x = 55, y = 52, width = 55, height=45,default.units = "mm")
plotGG(pirna,x = 105, y = 5,width = 75, height=110,default.units = "mm")

pageGuideHide()

dev.off()

