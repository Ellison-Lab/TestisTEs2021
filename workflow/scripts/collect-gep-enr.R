library(tidyverse)
library(arrow)

ont_df <- snakemake@input %>%
	  set_names(.,str_extract(.,'(?<=enrichment-).+(?=\\.csv)')) %>%
	  map_df(read_csv,.id='ont')

write_dataset(ont_df, snakemake@output[[1]], format = "arrow")
