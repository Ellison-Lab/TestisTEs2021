library(tidyverse)
library(jsonlite)
library(arrow)

clusters <- read_csv(snakemake@input[['clusters']]) %>%
  mutate_if(is.character,~str_remove_all(.,"[b\\']")) %>%
  mutate(clusters=as.character(clusters))

obsm <- read_csv(snakemake@input[['obsm']]) %>%
  select(contains('umap'))

obs_df <- read_csv(snakemake@input[['obs']]) %>%
  mutate_if(is.character,~str_remove_all(.,"[b\\']")) %>%
  left_join(clusters, by='clusters') %>%
  bind_cols(obsm)

write_dataset(obs_df, snakemake@output[[1]], format = "arrow")
