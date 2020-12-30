library(tidyverse)
library(jsonlite)
library(arrow)

message('opening var')
var_df <- open_dataset(snakemake@input[['var']], format='arrow') %>%
  dplyr::select(X1, gene_symbol) %>%
  collect()

message('opening obs')
obs_df <- open_dataset(snakemake@input[['obs']], format='arrow') %>%
    dplyr::select(X1, clusters, cell_type) %>%
    collect()

message('opening expr')
expr_df <- read_csv(snakemake@input[['expr']])

expr_df <- mutate(expr_df, index = str_extract(index,regex("(?<=b').+(?=')")))

message('gathering long')
expr_df <- expr_df %>% gather(gene_id, expression, -index)

message('adding gene symbol')
expr_df <- var_df %>% dplyr::select(X1,gene_symbol) %>%
    left_join(expr_df,., by=c(gene_id='X1'))

message('addition of cell type as partion')
# add cell types for partitioning
expr_df <- left_join(expr_df, obs_df, by=c(index="X1"))

message('partitioning')
# partion for quicker querying downstream
expr_df <- expr_df %>% group_by(clusters, cell_type)

write_dataset(expr_df, snakemake@output[[1]], format = "arrow")
