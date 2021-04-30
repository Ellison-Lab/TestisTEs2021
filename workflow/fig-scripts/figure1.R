library(tidyverse)

extrafont::loadfonts(quiet=TRUE)

library(patchwork)

source("workflow/fig-scripts/theme.R")

caption <- "Larval w1118 testis scRNA-seq cell type annotation. A) UMAP embedding of cells. Each cell is colored by its parent cluster, which corresponds to a cell type found in the larval testes. B) Dotplot shows average expression level (color) and the percentage of each cluster's cells expressing a given cell type marker."

umap <- read_rds('results/figs/intro_larval_umap/intro_larval_umap.ggp.rds') +
  theme(axis.text = element_text(size=rel(1)), axis.title = element_text(size=rel(1))) +
  theme(aspect.ratio = NULL)

markers <- read_rds('results/figs/larval_marker_expression/larval_marker_expression.ggp.rds') +
  theme(aspect.ratio = NULL) + 
  theme(axis.text = element_text(size=rel(0.5)),
        axis.text.x = element_text(size=rel(1)),
        axis.text.y = element_text(face="italic"),
        axis.title = element_text(size=rel(0.5)), 
        legend.key.size = unit(5,"pt"),legend.text = element_text(size=rel(0.3)), 
        legend.title = element_text(size=rel(0.3)))

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
  theme(axis.text.x = element_text(angle=90, hjust=1, size=rel(0.5)), axis.title.y =  element_text(size=rel(0.5)), axis.text.y=element_text(size=rel(0.5))) + xlab("") +
  ylab("cells per cluster") +
  guides(fill=F) + theme(axis.text = element_text(size=rel(1)), axis.title = element_text(size=rel(1))) +
  theme(aspect.ratio = NULL) +
  scale_y_continuous(breaks = seq(0,3000,by=500))

layout <-"
AAAABBB
AAAACCC
"

#cap <- "
#Figure 1: Identification of testis cell populations. A) UMAP projection groups transcriptionally similar cells in 2D space. Cells are colored by assigned cell type. B) Dot plot shows expression of selected marker genes used for cell type assignment. Color of each dot corresponds to mean normalized and log-transformed expression within cell clusters. Dot size corresponds to the proportion of cells in each cluster expressing the marker. C) Cell counts within each cell type cluster. Spermatocytes are the most abundant cell type by our analysis approach.
#"

p <- umap + markers + n_per_clust +
  plot_annotation(tag_levels = 'A',
                  #caption = str_wrap(cap,width = 200),
                  theme=theme(plot.caption = element_text(hjust=0, family="Arial"))) +
  plot_layout(design = layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5)))


ggsave(snakemake@output[[1]], p, width = 10, height = 6)

print(snakemake@output[[2]])

saveRDS(p,file=snakemake@output[[2]])
