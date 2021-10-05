library(tidyverse)
library(openxlsx)

source_data_list <- list(
f1a = read_tsv("results/figs/intro_larval_umap/intro_larval_umap.dat.tsv"),
f1b = read_tsv("results/figs/larval_marker_expression/larval_marker_expression.dat.tsv"),
f1c = read_tsv("results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.dat.tsv"),
f2a = read_tsv("results/figs/larval_te_heatmap/larval_te_heatmap.dat.tsv"),
f2b = read_tsv("results/figs/te_expression_by_cluster/te_expression_by_cluster.dat.tsv"),
f3a = read_tsv("results/figs/larval_all_gep_te_barchart/larval_all_gep_te_barchart.dat.tsv"),
f3b = read_tsv("results/figs/larval_tep_pie/larval_tep_pie.dat.tsv"),
f3c = read_tsv("results/figs/larval_tep_usage_umap/larval_tep_usage_umap.dat.tsv"),
f4a = read_tsv("results/figs/larval_fish_candidate_umis/larval_fish_candidate_umis.dat.tsv"),
f4c = read_tsv("results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.dat.tsv"),
f5a = read_tsv("results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.dat.tsv"),
f5b = read_tsv("results/figs/w1118_y_linked_copies/w1118_y_linked_copies.dat.tsv"),
f5c = read_tsv("results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.dat.tsv"),
f5d = read_tsv('results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.dat.tsv'),
f5e = read_tsv('results/figs/larval_pirna_expression/larval_pirna_expression.dat.tsv'),
st2 = read_tsv("results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.dat.tsv"), #rerun
s1a = tibble(), # rerun
s1b = read_tsv("results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.dat.tsv"),
s1c = read_tsv("results/figs/spermatogonia_marker_expr/spermatogonia_marker_expr.dat.tsv"),
s2a = read_tsv("results/figs/larval_raw_te_umis/larval_raw_te_umis.dat.tsv"),
s2b = read_tsv("results/figs/larval_normalized_te_umis/larval_normalized_te_umis.dat.tsv"),
s2c = read_tsv("results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.dat.tsv") %>%filter(str_detect(gene_id,"FBgn")),
s2d = read_tsv("results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.dat.tsv") %>% filter(!str_detect(gene_id,"FBgn")),
s2e = readRDS("results/figs/w1118_full_length_tx/w1118_full_length_tx.line.ggp.rds")$data,
s2f = readRDS("results/figs/w1118_full_length_tx/w1118_full_length_tx.heat.ggp.rds")$data,
s2g = readRDS("results/figs/chimerics/chimerics.ggp.rds")$data,
s2h = readRDS("results/figs/chimerics/chimerics.supporting_reads.ggp.rds")$data,
s3a = readRDS("results/figs/ica_reproducibility/ica_reproducibility.pval.ggp.rds")$data,
s3b = read_tsv("results/figs/larval_gep_size_histogram/larval_gep_size_histogram.dat.tsv"),
s3c = read_tsv("results/figs/larval_ica_optimized_enr/larval_ica_optimized_enr.dat.tsv"),
s7 = read_rds("results/figs/all_dataset_tep_scores/all_dataset_tep_scores.tes.ggp.rds")$data,
s8a = read_rds("results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.tpaf.ggp.rds")$data,
s8b = read_rds("results/figs/larval_tep_go_enrichment/larval_tep_go_enrichment.ggp.rds")$data,
s8c = read_tsv("results/figs/larval_y_gene_enr/larval_y_gene_enr.dat.tsv"), # rerun
s8d = read_tsv("results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.dat.tsv"),
s8e = read_rds("results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.box.paired.ggp.rds")$data,
s9a = read_tsv("results/figs/flam_tep_enrichment/flam_tep_enrichment.dat.tsv"), # reru
s9b = ("results/figs/ovary_silenced_enrichment/ovary_silenced_enrichment.dat.tsv"), #rerun
s10 = read_rds("results/figs/arriba_total_rna_fusions/arriba_total_rna_fusions.ggp.rds")$data,
s11 = read_rds("results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.4.ggp.rds")$data,
s12 = read_rds("results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.batch.ggp.rds")$data,
s13a = read_tsv("results/figs/larval_gep_corr_heatmap/larval_gep_corr.dat.tsv") %>% filter(calculated.from == "by.usage.score"),
s13b = read_tsv("results/figs/larval_gep_corr_heatmap/larval_gep_corr.dat.tsv") %>% filter(calculated.from == "by.membership.score"))



wb <- createWorkbook()

for (x in names(source_data_list)) {
  addWorksheet(wb,sheetName = x)  
  writeData(wb, sheet = x, x = source_data_list[[x]])
}

saveWorkbook(wb, "~/Downloads/2110051000_SourceData.xlsx", overwrite = TRUE)
