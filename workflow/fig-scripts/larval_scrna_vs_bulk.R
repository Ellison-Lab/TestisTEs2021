library(tidyverse)
library(ggpubr)
library(ggpointdensity)
library(arrow)
library(ragg)

source("workflow/fig-scripts/theme.R")

combine_transposon_bits <- function(tnm,mat) {
  raw_tn <- tnm[str_detect(tnm,"-I") | str_detect(tnm,"-LTR") | str_detect(tnm,"_I") | str_detect(tnm,"_LTR")]
  
  summed_mat <- raw_tn %>% 
    split(.,str_extract(.,regex("[[:alnum:]]+(?=[-_])"))) %>%
    map(.f=function(x) mat[x,]) %>%
    map(as.matrix) %>%
    map(.f=function(x) {
      as_tibble(x) %>% 
        gather(cell,count) %>% 
        group_by(cell) %>% 
        summarize(count=sum(count)) %>% 
        spread(cell,count) %>% as.matrix()
    }) %>%
    enframe() %>%
    mutate(nc = map_int(value,ncol)) %>% 
    filter(nc==ncol(mat)) %>%
    dplyr::select(name,value) %>%
    deframe() %>% {
      a <- .
      do.call(rbind,a) %>%
        magrittr::set_rownames(value = names(a))
    } %>%
    Matrix::Matrix(.,sparse=T)
  
  rbind(mat[!rownames(mat) %in% raw_tn,],summed_mat)
}

larv_bulk_01 <- read_tsv("results/finalized/larval-polya/larval_testes_cleaned_papain_01.tsv")
larv_bulk_02 <- read_tsv("results/finalized/larval-polya/larval_testes_papain_02.tsv")
larv_bulk_03 <- read_tsv("results/finalized/larval-polya/larval_testes_papain_03.tsv")
larv_bulk_04 <- read_tsv("results/finalized/larval-polya/larval_testes_papain_04.tsv")

larv_bulk <- bind_rows(list(`bulk rep 1`=larv_bulk_01, `bulk rep 2`=larv_bulk_02, `bulk rep 3`=larv_bulk_03, `bulk rep 4` =larv_bulk_04), .id="replicate") %>%
  gather(logic, count, -replicate, -gene_id)

larv_bulk <- filter(larv_bulk , 
                    !gene_id %in% c("N_unmapped","N_multimapping", "N_noFeature", "N_ambiguous")) %>%
  filter(logic == "second_strand")

trnrnas <- larv_bulk %>% filter(!str_detect(gene_id,"FBgn")) %>% pull(gene_id)

# counts - sum together LTR te portions
larv_bulk <- larv_bulk %>% dplyr::select(-logic) %>%
  spread(replicate,count) %>%
  column_to_rownames("gene_id") %>%
  as.matrix() %>%
  combine_transposon_bits(trnrnas,.) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

# counts
pseudobulk <- open_dataset("results/finalized/larval-w1118-testes/expr/", format='arrow') %>%
  collect() %>%
  mutate(scRNA = exp(expression) - 1) %>%
  group_by(gene_id) %>%
  summarize(scRNA = sum(scRNA))

df <- full_join(larv_bulk, pseudobulk, by="gene_id")

df <- df %>%
  mutate_if(is.numeric, ~replace_na(.,0))

df <- df %>%
  gather(bulk.replicate, bulk, -gene_id, -scRNA) %>%
  mutate(bulk = log2(bulk + 1)) %>%
  mutate(scRNA = log2(scRNA + 1))

g1 <- ggplot(df,aes(scRNA,bulk)) +
  facet_wrap(~bulk.replicate, scales="free") +
  geom_pointdensity()+
  stat_cor(label.sep = "\n", method='spearman') +
  #guides(color=F) +
  theme_gte21() +
  theme(aspect.ratio = 1) +
  scale_color_gte21("diverging", discrete = F, limits = c(0,500), oob=scales::squish, name="Density")

g2 <- df %>%
  filter(!str_detect(gene_id,"FBgn")) %>%
ggplot(aes(scRNA,bulk)) +
  facet_wrap(~bulk.replicate, scales="free") +
  geom_point(color='red') +
  stat_cor(label.sep = "\n", method='spearman') +
  #guides(color=F) +
  theme_gte21() +
  theme(aspect.ratio = 1) +
  scale_color_gte21("diverging", discrete = F, limits = c(0,500), oob=scales::squish, name="Density")


agg_png(snakemake@output[['png1']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g1)
dev.off()

agg_png(snakemake@output[['png2']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()

saveRDS(g1,snakemake@output[['ggp1']])
saveRDS(g2,snakemake@output[['ggp2']])

write_tsv(df,snakemake@output[['dat']])


# Export stats info -----------------------------------------------------------------------------------

raw.stats <- list(`all features` = g1, `only TEs`=g2) %>%
  map(~.$data) %>%
  map(~split(.,.$bulk.replicate)) %>%
  map(~map_df(.,~{broom::tidy(cor.test(x=.$scRNA, y=.$bulk,method = "spearman"),)},.id="comparison")) %>%
  bind_rows(.id="comparison0")


stats.export <- raw.stats %>%
  mutate(script= "larval_scrna_vs_bulk.R") %>%
  mutate(desc = "correlation of bulk replicates with scRNA") %>%
  mutate(func = "stats::cor.test/ggpubr::stat_cor") %>%
  mutate(ci = NA) %>%
  mutate(comparison = paste(comparison,"vs scRNA using",comparison0)) %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic=estimate, ci)

write_tsv(stats.export,snakemake@output[['stats']])
