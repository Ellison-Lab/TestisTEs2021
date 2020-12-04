library(rtracklayer)
library(tidyverse)

log <- file(snakemake@log[[1]], open="wt")
sink(log)

fa <- rtracklayer::import(snakemake@input[[1]])

nms <- names(fa) %>% str_extract(regex(".+(?=#)"))

nms %>% is.na %>% any %>% {!.} %>% stopifnot

names(fa) <- nms

rtracklayer::export(fa, snakemake@output[[1]])
