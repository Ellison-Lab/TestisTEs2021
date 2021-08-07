library(tidyverse)
library(jsonlite)
library(arrow)
library(ragg)

source("workflow/fig-scripts/theme.R")

te_lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")%>%
  dplyr::select(merged_te, te=gene_id) %>%
  distinct()

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership", format='arrow')

tep.name <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  left_join(te_lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  distinct() %>%
  group_by(module) %>%
  summarize(n_tes=n()) %>%
  arrange(-n_tes) %>%
  head(1) %>%
  pull(module) %>% as.character()

tep.members <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  filter(as.character(module)==tep.name) %>%
  pull(X1)


arriba_df <- Sys.glob("results/finalized/arriba/w1118_testes*fusions.tsv") %>%
  set_names(.,str_extract(.,"(?<=arriba\\/).+(?=\\.fusions.tsv)")) %>%
  map_df(read_tsv,.id = "replicate")

#arriba_df <- Sys.glob("subworkflows/gte21-total-rna-fusions/results/arriba/w1118_testes_*/fusions.tsv") %>%
#  set_names(.,str_extract(.,"(?<=arriba\\/).+(?=\\/)")) %>%
#  map_df(read_tsv,.id = "replicate")

colnames(arriba_df) <- arriba_df %>%
  names() %>%
  str_remove_all(.,"#") %>%
  tidy_names(syntactic = T, quiet = T)

arriba_df <- dplyr::select(arriba_df,-read_identifiers)


# Make a df with the parsed TE fusions. I identify and rename columns to be easier to interpret as relating to the TE
# or the gene. No filtering beyond each entry must have 1 gene, 1 TE.

te_gene_fusion_df <- arriba_df %>%
  filter(xor(str_detect(gene1,"FBgn"),str_detect(gene2,"FBgn"))) %>%
  filter(gene1 != "." & gene2 != ".") %>%
  mutate(te = ifelse(str_detect(gene1,"FBgn"),gene2,gene1)) %>%
  mutate(partner = ifelse(str_detect(gene1,"FBgn"),gene1,gene2)) %>%
  mutate(te.strand = ifelse(str_detect(gene1,"FBgn"),strand2.gene.fusion.,strand1.gene.fusion.)) %>%
  mutate(partner.strand = ifelse(str_detect(gene1,"FBgn"),strand1.gene.fusion.,strand2.gene.fusion.)) %>%
  mutate(te.breakpoint = ifelse(str_detect(gene1,"FBgn"),breakpoint2,breakpoint1)) %>%
  mutate(partner.breakpoint = ifelse(str_detect(gene1,"FBgn"),breakpoint1,breakpoint2)) %>%
  mutate(te.coverage = ifelse(str_detect(gene1,"FBgn"),coverage2,coverage1)) %>%
  mutate(partner.coverage = ifelse(str_detect(gene1,"FBgn"),coverage1,coverage2)) %>%
  mutate(te.split_reads = ifelse(str_detect(gene1,"FBgn"),split_reads2,split_reads1)) %>%
  mutate(partner.split_reads = ifelse(str_detect(gene1,"FBgn"),split_reads1,split_reads2)) %>%
  mutate(partner.site = ifelse(str_detect(gene1,"FBgn"),site1,site2)) %>%
  dplyr::select(replicate, te, partner, te.strand, partner.strand, te.breakpoint, partner.breakpoint, partner.site, type, tags, reading_frame, confidence, te.coverage, partner.coverage, te.split_reads, partner.split_reads, discordant_mates)


# Next map to the merged IDs and annotate by TEP membership.

stopifnot(nrow(te_gene_fusion_df[!te_gene_fusion_df$te %in% te_lookup$te,]) == 0 )

te_gene_fusion_df <-  mutate(te_gene_fusion_df, te.supporting_reads = te.split_reads + discordant_mates, partner.supporting_reads = partner.split_reads + discordant_mates)

te_gene_fusion_df <- filter(te_gene_fusion_df, partner.supporting_reads > 0 & te.supporting_reads > 0)

te_gene_fusion_df <- left_join(te_gene_fusion_df,te_lookup)

te_gene_fusion_df <- mutate(te_gene_fusion_df, te.gep = ifelse(merged_te %in% tep.members,"TEP","other")) %>%
  mutate(te_gene_fusion_df, partner.gep = ifelse(partner %in% tep.members,"TEP","other"))

rep_exact_support_df <- te_gene_fusion_df %>%
  dplyr::select(replicate,te.breakpoint, partner.breakpoint) %>%
  distinct() %>%
  group_by(te.breakpoint, partner.breakpoint) %>%
  tally(name = "n_reps_exact") %>% arrange(-n_reps_exact) %>%
  ungroup()

rep_partner_support_df <- te_gene_fusion_df %>%
  group_by(merged_te, partner) %>%
  summarize(n_reps_partners = length(unique(replicate))) %>% arrange(-n_reps_partners) %>%
  ungroup()

te_gene_fusion_df <- left_join(te_gene_fusion_df, rep_exact_support_df) %>% left_join(rep_partner_support_df)

te_gene_fusion_df <- te_gene_fusion_df %>%
  mutate(partner.chrom = str_extract(partner.breakpoint,".+(?=:)")) %>%
  mutate(te.percent.supporting = te.supporting_reads / (te.supporting_reads + te.coverage),
         partner.percent.supporting = partner.supporting_reads / (partner.supporting_reads + partner.coverage))

g <- te_gene_fusion_df %>%
  filter(n_reps_exact > 1) %>%
  group_by(te, partner) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  filter(!str_detect(partner,",")) %>%
  filter(partner.chrom=="Y") %>%
  filter(te.gep == "TEP") %>%
  mutate(id = paste(partner, merged_te,  sep=" // ")) %>%
  mutate(partner.gep = paste("fused to",partner.gep,"gene")) %>%
  ggplot(aes(reorder(str_replace_all(id, "_(?=[Y])","\n"), te.percent.supporting),te.percent.supporting,fill=replicate)) +
  geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("breakpoint") +
  ylab("% supporting") +scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) + theme_gte21() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  facet_grid(.~partner.gep, scales = "free_x" ) +
  scale_fill_gte21()


agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])

write_tsv(te_gene_fusion_df,snakemake@output[['dat']])
