library(LoomExperiment)
library(monocle3)
library(tidyverse)
library(org.Dm.eg.db)
library(garnett)
library(gt)

set.seed(snakemake@params[['seed']])

fl <- "wf/preproc-testis/out/larval-testes.ingest-integrated.raw-layer.loom"
fl <- snakemake@input[['larval']]

obs_fl <- "wf/preproc-testis/out/larval-testes.anno/obs.csv"
obs_fl <- paste0(snakemake@input[['csvs']],'/obs.csv')

var_fl <- "wf/preproc-testis/out/larval-testes.anno/var.csv"
var_fl <- paste0(snakemake@input[['csvs']],'/var.csv')

markers_fl <- snakemake@input[['markers']]
ofl <- snakemake@output[["ofl"]]
g_marker_check_ofl <- snakemake@output[["g_marker_check"]]
g_table_ofl <- snakemake@output[["g_table"]]
g_garnett_results_ofl <- snakemake@output[["g_garnett_results"]]

cores <- snakemake@threads[[1]]

obs <- read_csv(obs_fl) %>%
	mutate_if(is.character,~str_extract(.,regex("(?<=b').+(?=')"))) %>%
	column_to_rownames("X1")

var <- read_csv(var_fl) %>%
	mutate_if(is.character,~str_extract(.,regex("(?<=b').+(?=')"))) %>%
	mutate(gene_symbol = ifelse(is.na(gene_symbol),annotation_ID,gene_symbol)) %>%
  mutate(gene_short_name = `gene_symbol`) %>%
  column_to_rownames('gene_symbol')

scle <- import(fl, type="SingleCellLoomExperiment")

X <- assay(scle,'matrix')
colnames(X) <- rownames(obs)

rownames(X) <- rownames(var)

pd <- new("AnnotatedDataFrame", data = obs)
fd <- new("AnnotatedDataFrame", data = var)

cds <- newCellDataSet(as(X, "dgCMatrix"),
                             phenoData = pd,
                             featureData = fd)

cds <- estimateSizeFactors(cds)

marker_check <- check_markers(cds, markers_fl,
                              db=org.Dm.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

g_marker_check <- plot_markers(marker_check)

classifier <- train_cell_classifier(cds = cds,cores=cores,
                                         marker_file = markers_fl,
                                         db=org.Dm.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 500,
                                         marker_file_gene_id_type = "SYMBOL")

cds <- classify_cells(cds, classifier,
                           db = org.Dm.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "SYMBOL")

res <- cds@phenoData@data %>%
  as_tibble(rownames = 'X1') %>%
  dplyr::select(X1,clusters, garnett_cluster,cell_type, cluster_ext_type) %>%
  mutate_at(vars('clusters','garnett_cluster'), as.character)

g_garnett_results <- res %>%
  ggplot(aes(clusters, fill=cell_type)) +
  geom_bar(position='fill', width=1, color='black') +
  theme_minimal() +
  scale_fill_viridis_d(option = 5, direction=-1) +
  theme(aspect.ratio = 1) +
  ylab("Percentage") +
  xlab("Clusters")

# create table in the final form i will visualize
res2show <- res %>%
	group_by(clusters,cell_type) %>%
	summarize(n=n()) %>%
	ungroup() %>%
	mutate(clusters = paste0(" ",clusters)) %>%
	spread(clusters,n) %>%
	mutate_at(vars(-cell_type),replace_na,0)

# create initial gt object
g_table <-  res2show %>%
  gt(rowname_col = 'cell_type') %>%
  tab_header("Cell type assignments",
             subtitle = md("Predicted cell type counts in each Leiden community")) %>%
  tab_stubhead(md("Cell type")) %>%
  tab_spanner(columns = everything(), label=md("`Leiden`"))

# color cells in gt object by top 2 represented cell types per
# community
for (l in paste0(" ",unique(res$clusters))) {
  print(sym(l))
  g_table <- g_table %>%
    #tab_style(style=cell_fill(color="gold"), locations = cells_body(columns = vars(l)))
    tab_style(style=cell_fill(color="goldenrod"),
              locations = cells_body(columns=vars(l),
                                     rows= res2show[[l]] %in% sort(res2show[[l]],T)[1:2]))
}

# color the unknown row light gray
g_table <- g_table %>%
  tab_style(style=cell_fill(color="lightgray"),
            locations = cells_body(columns=everything(), rows=cell_type == "Unknown"))

# output final summary table
res %>%
  group_by(clusters,cell_type) %>%
  summarize(n=n()) %>%
  group_by(clusters) %>%
  mutate(pct = n/sum(n)) %>%
  filter(cell_type != "Unknown") %>%
  top_n(1,n) %>%
  ungroup() %>%
  mutate(clusters2 = paste(clusters,cell_type,sep="/")) %>%
  write_csv(ofl)

# output plots
ggsave(g_marker_check_ofl, g_marker_check)
gtsave(g_table, g_table_ofl)
ggsave(g_garnett_results_ofl, g_garnett_results)
