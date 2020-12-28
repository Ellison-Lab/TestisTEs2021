library(tidyverse)
library(jsonlite)
library(arrow)

var_df <- read_csv(snakemake@input[['var']]) %>%
  mutate_if(is.character,~str_remove_all(.,"[b\\']"))

expr_df <- read_csv(snakemake@input[['expr']]) %>%
  mutate_if(is.character,~str_remove_all(.,"[b\\']"))

expr_df <- expr_df %>% gather(gene_id, expression, -index)

expr_df <- var_df %>% dplyr::select(X1,gene_symbol) %>%
    left_join(expr_df,., by=c(gene_id='X1'))

write_dataset(expr_df, snakemake@output[[1]], format = "arrow")
