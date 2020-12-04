library(rtracklayer)
library(tidyverse)

log <- file(snakemake@log[[1]], open="wt")
sink(log)

gr <- rtracklayer::import(snakemake@input[[1]])

# cellranger gets angry when strand != {+,-}
strand(gr[strand(gr) == "*"]) <- "+"

# make into a tbl for easy manipulation
# also fill in transcript symbol and enough fields to make cellranger happy
# and comfortable making an index from this
# lastly, feature type must be exon to be counted
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
gr <- gr %>%
  as_tibble() %>%
  mutate(seqnames=as.character(seqnames)) %>%
  mutate(gene_id = ifelse(is.na(gene_id),seqnames,gene_id)) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol),seqnames,gene_symbol)) %>%
  mutate(transcript_symbol = ifelse(is.na(transcript_symbol),seqnames,transcript_symbol)) %>%
  mutate(transcript_id = ifelse(is.na(transcript_id),seqnames,transcript_id)) %>%
  mutate(type = as.character(type)) %>%
  mutate(type = ifelse(type=="sequence_feature","exon",type))

# get the rna transcript types we don't want from the snakmake object
disallow <- gr %>%
   filter(type %in% snakemake@params[["disallow"]]) %>%
   dplyr::select(transcript_symbol, type)

# remove any lines corresponding to txs in these categories
gr <- gr %>% filter(!(transcript_symbol %in% disallow$transcript_symbol))

# convert to gr for export
gr <- GRanges(gr)

rtracklayer::export(gr, snakemake@output[[1]])
