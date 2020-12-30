library(tidyverse)
library(jsonlite)
library(arrow)

var_df <- read_csv(snakemake@input[[1]]) %>%
	mutate_if(is.character,~str_extract(.,regex("(?<=b').+(?=')")))

write_dataset(var_df, snakemake@output[[1]], format = "arrow")
