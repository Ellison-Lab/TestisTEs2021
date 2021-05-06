library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(readxl)
library(arrow)
library(ragg)

source("workflow/fig-scripts/theme.R")
chimerics <- Sys.glob("~/amarel-scratch/TE-proj-reorg/gte21-chimeric-rnaseq/results/breakpoints/*.breakpoint-depths.csv") %>%
  map_df(read_csv)

lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")

gtf <- import("~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-custom-genome/results/custom-genome/combined.fixed.gtf")

genes <- gtf[gtf$type == "mRNA"]

main_chroms <- c("2L", "2R","3L","3R","4","X","Y")

# to approximate cellranger counting, only count sense tx. This library has '-' sense reads.
chimerics <- chimerics %>% filter(other_strand == "-") 

res <- chimerics %>%
  filter(other_chr %in% lookup$gene_id) %>%
  mutate(end = start) %>% GRanges() %>%
  subsetByOverlaps(genes, ignore.strand=T) %>%
  as_tibble()

to_plot <- res %>% 
  mutate(pct.supporting = supporting_reads/(supporting_reads + other_reads)) %>%
  distinct() %>%
  group_by(breakpoint_id) %>%
  filter(length(unique(sample)) >=2) %>%
  filter(all(supporting_reads > 1)) %>%
  arrange(breakpoint_id) %>%
  ungroup()

g <- to_plot %>% dplyr::select(breakpoint_id, other_chr) %>%
  distinct() %>%
  left_join(distinct(lookup %>% dplyr::select(merged_te, gene_id)), by=c(other_chr="gene_id")) %>%
  dplyr::select(breakpoint_id, merged_te) %>%
  distinct() %>%
  group_by(merged_te) %>%
  tally() %>%
  ggplot(aes(reorder(merged_te,n), n)) +
  geom_col() +
  theme_gte21() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) +
  xlab("") + ylab("reproducible genic breakpoints") +
  scale_y_continuous(expand = expand_scale(mult = c(0,0), add = c(0,1)))

g2 <- to_plot %>%
  dplyr::select(breakpoint_id, other_chr,pct.supporting, sample, supporting_reads) %>%
  distinct() %>%
  left_join(distinct(lookup %>% dplyr::select(merged_te, gene_id)), by=c(other_chr="gene_id")) %>%
  dplyr::select(breakpoint_id, merged_te, sample,supporting_reads, pct.supporting) %>%
  ggplot(aes(reorder(merged_te,supporting_reads),supporting_reads)) +
  #stat_summary(position = "dodge") +
  geom_jitter(width = 0.1) +
  facet_wrap(~sample) +
  theme_gte21() +
  xlab("") + ylab("supporting reads") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

agg_png(snakemake@output[['png2']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
saveRDS(g2,snakemake@output[['ggp2']])


