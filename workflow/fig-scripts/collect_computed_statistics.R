library(tidyverse)

#fls <- Sys.glob("results/figs/*/*stats.tsv")
fls <- snakemake@input

df <- fls %>%
  map_df(read_tsv)
  
write_tsv(df,snakemake@output[["stats"]])
