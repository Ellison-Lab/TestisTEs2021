library(tidyverse)
library(readxl)
library(arrow)
library(ragg)

source("workflow/fig-scripts/theme.R")

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')
w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/expr", format='arrow')

df <- readxl::read_xlsx("resources/41467_2021_20897_MOESM5_ESM.xlsx",sheet = 'Gene Level Data') %>%
  dplyr::select(FBgn, contains('tpm')) %>%
  drop_na()

df.cells <- readxl::read_xlsx("resources/41467_2021_20897_MOESM5_ESM.xlsx",sheet = 'Cell Level Data')

# collect our data
our.data <- w1118.expr %>% 
  filter(gene_id %in% df$FBgn) %>% 
  collect() %>%
  unite(clusters2, clusters,cell_type,sep='/') %>%
  left_join(rename.table)

our.data.tpm.clusters <- our.data %>% 
  mutate(tpm = exp(expression) - 1) %>%
  #left_join(collect(w1118.obs), by=c(index='X1'))
  group_by(clusters.rename, gene_id) %>%
  summarise(tpm = mean(tpm), .groups = 'drop')

df2 <- df %>% gather(cluster.mahadevajaru, tpm,-FBgn) %>%
  separate(cluster.mahadevajaru, into = c('metric','clusters.mahadevajaru','rep'),extra = 'merge') %>%
  group_by(FBgn, clusters.mahadevajaru) %>%
  summarize(tpm = sum(tpm), .groups = 'drop') %>%
  dplyr::select(gene_id=FBgn, clusters.mahadevajaru, tpm)

cor.df <- crossing(this.study=unique(our.data.tpm.clusters$clusters.rename), mahadevajaru =unique(df2$clusters.mahadevajaru)) %>%
  mutate_all(as.character) %>%
  mutate(clusters = map2(.x=this.study, .y=mahadevajaru, ~list(.x,.y))) %>%
  pull(clusters) %>%
  set_names(.,map_chr(.,~paste(.[1],.[2], sep="_"))) %>%
  map(.f=~inner_join(filter(our.data.tpm.clusters,clusters.rename == .x[1]),filter(df2, clusters.mahadevajaru == .x[2]),by='gene_id')) %>%
  map_df(~broom::tidy(cor.test(.$tpm.x,.$tpm.y, method='spearman')), .id='pair') %>%
  separate(pair, into=c("x",'y'), sep="_")

g <- cor.df %>%
  dplyr::select(x,y, estimate) %>%
  mutate(x=fct_reorder(x,as.numeric(str_extract(x,"\\d+")))) %>%
  mutate(y=fct_relevel(y,c('G','E1','M1','L1','C1','C2','C3','C4','P','T'))) %>%
  #spread(y, estimate) %>%
  ggplot(aes(y,x,fill=estimate)) +
  geom_tile() +
  scale_fill_distiller(type = 'seq', palette = 8, direction = 1,name='Spearman') +
  theme_minimal() +
  coord_fixed() +
  theme(aspect.ratio = 1) +
  xlab("Mahadevaraju et al. clusters") + ylab("This study")


n_mahadevaraju <- df.cells %>% 
  mutate(is_cell_used_in_study = as.logical(is_cell_used_in_study)) %>%
  filter(is_cell_used_in_study) %>%
  tally()

n_this_study <- w1118.obs %>% collect() %>% tally()

n_cells <- bind_rows(`Mahadevaraju et al.` = n_mahadevaraju, `This study` = n_this_study, .id = "study")


g2 <- ggplot(n_cells,aes(study, n)) +
  geom_col(aes(fill=study)) +
  theme_gte21() +
  scale_fill_gte21("binary", reverse = T) +
  guides(fill=F) + xlab("")


agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

agg_png(snakemake@output[['png2']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(2)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
saveRDS(g2,snakemake@output[['ggp2']])
write_tsv(cor.df,snakemake@output[['dat']])

