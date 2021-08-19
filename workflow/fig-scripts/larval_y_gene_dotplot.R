library(tidyverse)
library(arrow)
library(ragg)
library(rtracklayer)

source("workflow/fig-scripts/theme.R")

gtf <- import('subworkflows/gte21-custom-genome/results/custom-genome/combined.fixed.gtf') %>%
  as_tibble()

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')
w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/expr/", format='arrow')

ylinked <- unique(gtf %>% filter(seqnames == 'Y' & type == 'mRNA') %>% pull(gene_symbol)) %>%
  tibble(gene_symbol = ., group='y-linked')

tmac <- c('aly','wuc','tomb',"kmg") %>% tibble(gene_symbol = ., group='tMAC')
ttaf <- c('sa','nht','mia','can','Taf12L') %>% tibble(gene_symbol = ., group='tTAF')
tbrd <- c('tbrd-1','tbrd-2') %>% tibble(gene_symbol = ., group='tBRD')
tplus <- c('tplus3a','tplus3b',"bol") %>% tibble(gene_symbol = ., group='tPAF')
tep.marker <- c('EAChm') %>% tibble(gene_symbol = ., group='TEP-HI')

male.meiosis1.associated <- bind_rows(ylinked, tmac, ttaf, tbrd, tplus, tep.marker)

markers2display <- male.meiosis1.associated$gene_symbol

dat <- map_df(w1118.obs %>% collect() %>% pull(clusters) %>% unique() %>% as.list %>% set_names(.,.),
              ~{filter(w1118.expr, clusters == . & gene_symbol %in% markers2display) %>% collect()}) %>%
  left_join(collect(w1118.obs), by=c(index='X1','cell_type'='cell_type')) %>%
  left_join(rename.table) %>%
  ungroup() %>%
  #group_by(cell_type) %>%
  mutate(ord = dense_rank(gene_symbol)) %>%
  mutate(gene_symbol = fct_relevel(gene_symbol, markers2display))

dat2 <- dat %>%
  group_by(clusters.rename, gene_symbol) %>%
  summarize(pct.expressing = sum(expression > 0)/n(), mean.expression=mean(expression)) %>%
  left_join(male.meiosis1.associated)

g <- dat2 %>%
  filter(group == "y-linked") %>%
  ggplot(aes(gene_symbol, clusters.rename)) +
  geom_point(aes(size=pct.expressing, fill=mean.expression), shape=21) +
  scale_fill_fermenter(palette = 8, direction = 1, name='mean(log-norm UMIs)', guide=guide_legend(label.position = 'bottom', title.position = 'top')) +
  scale_size(range=c(0, rel(7)), name='Proportion expressing', guide=guide_legend(label.position = 'bottom', title.position = 'top')) +
  theme_gte21()  +
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.text.y=element_text(face="italic")) +
  theme(legend.text = element_text(size=rel(0.5)), legend.title = element_text(size=rel(0.5))) +
  theme(aspect.ratio = 1.5) +
  coord_flip() +
  xlab('') + ylab('') +
  facet_wrap(~group, scales = "free", ncol = 1)


g2 <-  dat2 %>%
  filter(group %in% c("tBRD","tPAF","tTAF","tMAC")) %>%
  ggplot(aes(gene_symbol, clusters.rename)) +
  geom_point(aes(size=pct.expressing, fill=mean.expression), shape=21) +
  scale_fill_fermenter(palette = 8, direction = 1, name='mean(log-norm UMIs)', guide=guide_legend(label.position = 'bottom', title.position = 'top')) +
  scale_size(range=c(0, rel(3.2)), name='Proportion expressing', guide=guide_legend(label.position = 'bottom', title.position = 'top')) +
  theme_gte21()  +
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.text.y=element_text(face="italic")) +
  theme(legend.text = element_text(size=rel(0.5)), legend.title = element_text(size=rel(0.5))) +
  theme(aspect.ratio = 2) +
  coord_flip() +
  xlab('') + ylab('')

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

agg_png(snakemake@output[['png2']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()


saveRDS(g,snakemake@output[['ggp']])
saveRDS(g2,snakemake@output[['ggp_tpaf']])
write_tsv(dat2,snakemake@output[['dat']])
