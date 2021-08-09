library(tidyverse)

#celltypes <- read_tsv("results/figs/celltype_rename_table.tsv")
celltypes <- read_tsv(snakemake@input[['rename']])

#df <- read_csv("../gte21-scrna/results/scanpy/larval-w1118-testes/diffs.csv.gz")
df <- read_csv(snakemake@input[['diffs']])

res <- celltypes %>%
  mutate(id = as.numeric(str_extract(clusters2,"\\d+(?=\\/)"))) %>%
  left_join(df) %>%
  dplyr::select(clusters.orig = clusters2, clusters.rename, clusters_id.orig = id, gene_id = names, logfoldchanges, pvals, pvals_adj)

write_tsv(res, snakemake@output[["tsv"]])
