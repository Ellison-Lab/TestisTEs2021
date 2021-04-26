library(tidyverse)
library(readxl)
library(arrow)
library(ragg)

source("workflow/fig-scripts/theme.R")

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')

obs <- w1118.obs %>% collect() %>%
  left_join(rename.table)

obs <- obs %>% dplyr::select(X1, batch, doublet_score, n_genes, log1p_total_counts, percent_mito, phase, clusters.rename)

g.dubs <- obs %>%
  ggplot(aes(clusters.rename, doublet_score, fill=clusters.rename)) +
  geom_violin() +
  theme_gte21() +
  theme(aspect.ratio = 0.5, axis.text.x = element_text(angle=90, hjust=1)) +
  scale_fill_gte21() +
  guides(fill=F) +
  xlab("") + ylab('Scrublet score (post-filtering)')

g.batch <- obs %>%
  ggplot(aes(clusters.rename, fill=batch)) +
  geom_bar(position='stack') +
  theme_gte21() +
  theme(aspect.ratio = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_brewer(type='qual', palette = 3) +
  xlab("") + ylab('N cells')

g.n_genes <- obs %>%
  ggplot(aes(clusters.rename, n_genes, fill=clusters.rename)) +
  geom_violin() +
  theme_gte21() +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_fill_gte21() +
  guides(fill=F) +
  xlab("") + ylab('detected genes')

g.log1p_umis <- obs %>%
  ggplot(aes(clusters.rename, log1p_total_counts, fill=clusters.rename)) +
  geom_violin() +
  theme_gte21() +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_fill_gte21() +
  guides(fill=F) +
  xlab("") + ylab('log1p(UMIs)')

g.percent_mito <- obs %>%
  ggplot(aes(clusters.rename, percent_mito, fill=clusters.rename)) +
  geom_violin() +
  theme_gte21() +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_fill_gte21() +
  scale_y_continuous(labels=scales::percent) +
  guides(fill=F) +
  xlab("") + ylab('% Mitochondrial gene UMIs (post-filtering)')

g.phase <- obs %>%
  ggplot(aes(clusters.rename, fill=phase)) +
  geom_bar(position='stack') +
  theme_gte21() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_brewer(type='qual', palette = 3) +
  xlab("") + ylab('N cells')

agg_png(snakemake@output[['png_dubs']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g.dubs)
dev.off()

agg_png(snakemake@output[['png_batch']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g.batch)
dev.off()

agg_png(snakemake@output[['png_n_genes']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g.n_genes)
dev.off()

agg_png(snakemake@output[['png_log1p_umis']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g.log1p_umis)
dev.off()

agg_png(snakemake@output[['png_percent_mito']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g.percent_mito)
dev.off()

agg_png(snakemake@output[['png_phase']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g.phase)
dev.off()

saveRDS(g.dubs,snakemake@output[['ggp_dubs']])
saveRDS(g.batch,snakemake@output[['ggp_batch']])
saveRDS(g.n_genes,snakemake@output[['ggp_n_genes']])
saveRDS(g.log1p_umis,snakemake@output[['ggp_log1p_umis']])
saveRDS(g.percent_mito,snakemake@output[['ggp_percent_mito']])
saveRDS(g.phase,snakemake@output[['ggp_phase']])

write_tsv(obs,snakemake@output[['dat']])
