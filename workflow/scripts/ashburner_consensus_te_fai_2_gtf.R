library(rtracklayer)
library(tidyverse)

log <- file(snakemake@log[[1]], open="wt")
sink(log)

fai_tbl <- read_tsv(snakemake@input[[1]],
         col_names = c("name","length","offset","a","b"))


fai_tbl %>%
  mutate(start=1, end = length, seqnames=name) %>%
  dplyr::select(seqnames, start, end) %>%
  GRanges() %>%
  export(snakemake@output[[1]])
