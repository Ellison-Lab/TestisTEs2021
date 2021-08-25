rule celltype_rename_table:
    output:
        tsv="results/figs/celltype_rename_table.tsv"
    script:
        "../fig-scripts/celltype_rename_table.R"

rule intro_larval_umap:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/intro_larval_umap/intro_larval_umap.png",
        ggp = "results/figs/intro_larval_umap/intro_larval_umap.ggp.rds",
        dat = "results/figs/intro_larval_umap/intro_larval_umap.dat.tsv"
    script:
        "../fig-scripts/intro_larval_umap.R"

rule larval_marker_expression:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_marker_expression/larval_marker_expression.png",
        ggp = "results/figs/larval_marker_expression/larval_marker_expression.ggp.rds",
        dat = "results/figs/larval_marker_expression/larval_marker_expression.dat.tsv"
    script:
        "../fig-scripts/larval_marker_expression.R"


rule larval_scrna_vs_bulk:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png1 = "results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.all.png",
        png2 = "results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.tes.png",
        ggp1 = "results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.all.ggp.rds",
        ggp2 = "results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.tes.ggp.rds",
        dat = "results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.dat.tsv",
        stats = "results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.stats.tsv",
    script:
        "../fig-scripts/larval_scrna_vs_bulk.R"

rule larval_te_heatmap:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_te_heatmap/larval_te_heatmap.png",
        dat = "results/figs/larval_te_heatmap/larval_te_heatmap.dat.tsv"
    script:
        "../fig-scripts/larval_te_heatmap.R"


rule larval_normalized_te_umis:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_normalized_te_umis/larval_normalized_te_umis.png",
        ggp = "results/figs/larval_normalized_te_umis/larval_normalized_te_umis.ggp.rds",
        dat = "results/figs/larval_normalized_te_umis/larval_normalized_te_umis.dat.tsv",
        stats = "results/figs/larval_normalized_te_umis/larval_normalized_te_umis.stats.tsv"
    script:
        "../fig-scripts/larval_normalized_te_umis.R"

rule larval_raw_te_umis:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_raw_te_umis/larval_raw_te_umis.png",
        ggp = "results/figs/larval_raw_te_umis/larval_raw_te_umis.ggp.rds",
        dat = "results/figs/larval_raw_te_umis/larval_raw_te_umis.dat.tsv",
        stats = "results/figs/larval_raw_te_umis/larval_raw_te_umis.stats.tsv"
    script:
        "../fig-scripts/larval_raw_te_umis.R"

rule larval_later_sperm_marker_umis:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_later_sperm_marker_umis/larval_later_sperm_marker_umis.png",
        ggp = "results/figs/larval_later_sperm_marker_umis/larval_later_sperm_marker_umis.ggp.rds",
        dat = "results/figs/larval_later_sperm_marker_umis/larval_later_sperm_marker_umis.dat.tsv"
    script:
        "../fig-scripts/larval_later_sperm_marker_umis.R"


rule larval_tep_usage_umap:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_tep_usage_umap/larval_tep_usage_umap.png",
        ggp = "results/figs/larval_tep_usage_umap/larval_tep_usage_umap.ggp.rds",
        dat = "results/figs/larval_tep_usage_umap/larval_tep_usage_umap.dat.tsv"
    script:
        "../fig-scripts/larval_tep_usage_umap.R"

rule larval_all_gep_te_bachart:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_all_gep_te_barchart/larval_all_gep_te_barchart.png",
        ggp = "results/figs/larval_all_gep_te_barchart/larval_all_gep_te_barchart.ggp.rds",
        dat = "results/figs/larval_all_gep_te_barchart/larval_all_gep_te_barchart.dat.tsv"
    script:
        "../fig-scripts/larval_all_gep_te_barchart.R"


rule larval_tep_te_barchart:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_tep_te_barchart/larval_tep_te_barchart.png",
        ggp = "results/figs/larval_tep_te_barchart/larval_tep_te_barchart.ggp.rds",
        dat = "results/figs/larval_tep_te_barchart/larval_tep_te_barchart.dat.tsv"
    script:
        "../fig-scripts/larval_tep_te_barchart.R"

rule larval_tep_pie:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_tep_pie/larval_tep_pie.png",
        ggp = "results/figs/larval_tep_pie/larval_tep_pie.ggp.rds",
        dat = "results/figs/larval_tep_pie/larval_tep_pie.dat.tsv"
    script:
        "../fig-scripts/larval_tep_pie.R"

# --------------------------------------------


rule larval_cica_grid_heatmap:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_cica_grid_heatmap/larval_cica_grid_heatmap.png",
        ggp = "results/figs/larval_cica_grid_heatmap/larval_cica_grid_heatmap.ggp.rds",
        dat = "results/figs/larval_cica_grid_heatmap/larval_cica_grid_heatmap.dat.tsv"
    script:
        "../fig-scripts/larval_cica_grid_heatmap.R"

rule larval_gep_size_histogram:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_gep_size_histogram/larval_gep_size_histogram.png",
        ggp = "results/figs/larval_gep_size_histogram/larval_gep_size_histogram.ggp.rds",
        dat = "results/figs/larval_gep_size_histogram/larval_gep_size_histogram.dat.tsv"
    script:
        "../fig-scripts/larval_gep_size_histogram.R"

rule larval_tep_go_enrichment:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_tep_go_enrichment/larval_tep_go_enrichment.png",
        ggp = "results/figs/larval_tep_go_enrichment/larval_tep_go_enrichment.ggp.rds",
        dat = "results/figs/larval_tep_go_enrichment/larval_tep_go_enrichment.dat.tsv"
    script:
        "../fig-scripts/larval_tep_go_enrichment.R"

rule larval_ica_optimized_enr:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_ica_optimized_enr/larval_ica_optimized_enr.png",
        ggp = "results/figs/larval_ica_optimized_enr/larval_ica_optimized_enr.ggp.rds",
        dat = "results/figs/larval_ica_optimized_enr/larval_ica_optimized_enr.dat.tsv"
    script:
        "../fig-scripts/larval_ica_optimized_enr.R"

rule larval_gep_corr_heatmap:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png1 = "results/figs/larval_gep_corr_heatmap/larval_gep_corr.usage.png",
        png2 = "results/figs/larval_gep_corr_heatmap/larval_gep_corr.membership.png",
    script:
        "../fig-scripts/larval_gep_corr_heatmap.R"

rule larval_y_gene_enr:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_y_gene_enr/larval_y_gene_enr.png",
        ggp = "results/figs/larval_y_gene_enr/larval_y_gene_enr.ggp.rds",
        dat = "results/figs/larval_y_gene_enr/larval_y_gene_enr.dat.tsv",
        stats = "results/figs/larval_y_gene_enr/larval_y_gene_enr.stats.tsv"
    script:
        "../fig-scripts/larval_y_gene_enr.R"

rule larval_y_gene_dotplot:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.png",
        png2 = "results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.tpaf.png",
        ggp = "results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.ggp.rds",
        ggp_tpaf = "results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.tpaf.ggp.rds",
        dat = "results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.dat.tsv"
    script:
        "../fig-scripts/larval_y_gene_dotplot.R"

rule larval_fish_candidate_umis:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_fish_candidate_umis/larval_fish_candidate_umis.png",
        ggp = "results/figs/larval_fish_candidate_umis/larval_fish_candidate_umis.ggp.rds",
        dat = "results/figs/larval_fish_candidate_umis/larval_fish_candidate_umis.dat.tsv"
    script:
        "../fig-scripts/larval_fish_candidate_umis.R"

rule tidal_larval_tep_xa_boxplot:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png1 = "results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.png",
        ggp1 = "results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.ggp.rds",
        png2 = "results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.4.png",
        ggp2 = "results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.4.ggp.rds",
        dat = "results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.dat.tsv",
        stats = "results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.stats.tsv"
    script:
        "../fig-scripts/tidal-fig-v2.R"

rule larracuente_y_ins_barchart:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.png",
        ggp = "results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.ggp.rds",
        png2 = "results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.at_least_1.png",
        ggp2 = "results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.at_least_1.ggp.rds",
        dat = "results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.dat.tsv",
        stats = "results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.stats.tsv"
    script:
        "../fig-scripts/larracuente_y_ins_barchart.R"

rule larracuente_y_intronic_ins:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png1 = "results/figs/larracuente_y_intronic_ins/larracuente_y_intronic_ins.y_ins.png",
        png2 = "results/figs/larracuente_y_intronic_ins/larracuente_y_intronic_ins.y_ins_by_consensus.png",
        png3 = "results/figs/larracuente_y_intronic_ins/larracuente_y_intronic_ins.y_tep_ins.png",
        stats = "results/figs/larracuente_y_intronic_ins/larracuente_y_intronic_ins.y_tep_ins.stats.tsv",
    script:
        "../fig-scripts/larracuente_y_intronic_ins.R"


rule w1118_y_linked_copies:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png1 = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.1.png",
        ggp1 = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.1.ggp.rds",
        png2 = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.2.png",
        ggp2 = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.2.ggp.rds",
        png3 = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.3.png",
        ggp3 = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.3.ggp.rds",
        dat = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.dat.tsv",
        stats = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.stats.tsv"
    script:
        "../fig-scripts/w1118_y_linked_copies.R"

rule w1118_pct_y_linked_rna_vs_wgs:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.png",
        png2 = "results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.box.paired.png",
        ggp = "results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.ggp.rds",
        ggp2 = "results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.box.paired.ggp.rds",
        dat = "results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.dat.tsv",
        stats = "results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.stats.tsv"
    script:
        "../fig-scripts/w1118_pct_y_linked_rna_vs_wgs.R"

rule larval_pirna_genes:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_pirna_expression/larval_pirna_expression.png",
        ggp = "results/figs/larval_pirna_expression/larval_pirna_expression.ggp.rds",
        dat = "results/figs/larval_pirna_expression/larval_pirna_expression.dat.tsv"
    script:
        "../fig-scripts/larval_pirna_genes.R"

# ------------------------------------------------------------------------------

rule w1118_full_length_tx:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png1 = "results/figs/w1118_full_length_tx/w1118_full_length_tx.line.png",
        ggp1 = "results/figs/w1118_full_length_tx/w1118_full_length_tx.line.ggp.rds",
        png2 = "results/figs/w1118_full_length_tx/w1118_full_length_tx.heat.png",
        ggp2 = "results/figs/w1118_full_length_tx/w1118_full_length_tx.heat.ggp.rds",
        dat = "results/figs/w1118_full_length_tx/w1118_full_length_tx.dat.tsv",
        stats = "results/figs/w1118_full_length_tx/w1118_full_length_tx.stats.tsv"
    script:
        "../fig-scripts/w1118_full_length_tx.R"


rule comparison_with_mahadevaraju:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output,
        "resources/41467_2021_20897_MOESM5_ESM.xlsx"
    output:
        png = "results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.png",
        ggp = "results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.ggp.rds",
        png2 = "results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.n_cells.png",
        ggp2 = "results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.n_cells.ggp.rds",
        dat = "results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.dat.tsv",
        stats = "results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.stats.tsv"
    script:
        "../fig-scripts/comparison_with_mahadevaraju.R"

rule larval_scrna_basic_qc_stats:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output,
    output:
        png_dubs = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.dubs.png",
        ggp_dubs = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.dubs.ggp.rds",
        png_batch = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.batch.png",
        ggp_batch = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.batch.ggp.rds",
        png_n_genes = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.n_genes.png",
        ggp_n_genes = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.n_genes.ggp.rds",
        png_log1p_umis = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.log1p_umis.png",
        ggp_log1p_umis = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.log1p_umis.ggp.rds",
        png_percent_mito = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.percent_mito.png",
        ggp_percent_mito = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.percent_mito.ggp.rds",
        png_phase = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.phase.png",
        ggp_phase = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.phase.ggp.rds",
        dat = "results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.dat.tsv"
    script:
        "../fig-scripts/larval_scrna_basic_qc_stats.R"

rule all_dataset_tep_scores:
    input:
        rules.collect_mod_scores.output
    output:
        png_tes = "results/figs/all_dataset_tep_scores/all_dataset_tep_scores.tes.png",
        ggp_tes = "results/figs/all_dataset_tep_scores/all_dataset_tep_scores.tes.ggp.rds",
        dat_tes = "results/figs/all_dataset_tep_scores/all_dataset_tep_scores.tes.dat.tsv",
        stats = "results/figs/all_dataset_tep_scores/all_dataset_tep_scores.stats.tsv",
    script:
        "../fig-scripts/all_dataset_tep_scores.R"

rule chimerics:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output,
    output:
        png = "results/figs/chimerics/chimerics.png",
        ggp = "results/figs/chimerics/chimerics.ggp.rds",
        png2 = "results/figs/chimerics/chimerics.supporting_reads.png",
        ggp2 = "results/figs/chimerics/chimerics.supporting_reads.ggp.rds",
    script:
        "../fig-scripts/chimerics.R"

rule arriba_total_rna_fusions:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output,
    output:
        png = "results/figs/arriba_total_rna_fusions/arriba_total_rna_fusions.png",
        ggp = "results/figs/arriba_total_rna_fusions/arriba_total_rna_fusions.ggp.rds",
        dat = "results/figs/arriba_total_rna_fusions/arriba_total_rna_fusions.dat.tsv",
    script:
        "../fig-scripts/arriba_total_rna_fusions.R"

# ------------------------------------------------------------------------------
# Stuff for revisions
# ------------------------------------------------------------------------------

rule flam_tep_enrichment:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output,
    output:
        png = "results/figs/flam_tep_enrichment/flam_tep_enrichment.png",
        ggp = "results/figs/flam_tep_enrichment/flam_tep_enrichment.ggp.rds",
        dat = "results/figs/flam_tep_enrichment/flam_tep_enrichment.dat.tsv",
        stats = "results/figs/flam_tep_enrichment/flam_tep_enrichment.stats.tsv",
    script:
        "../fig-scripts/flam-enrich.R"

rule ovary_silenced_enrichment:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output,
    output:
        png = "results/figs/ovary_silenced_enrichment/ovary_silenced_enrichment.png",
        ggp = "results/figs/ovary_silenced_enrichment/ovary_silenced_enrichment.ggp.rds",
        dat = "results/figs/ovary_silenced_enrichment/ovary_silenced_enrichment.dat.tsv",
        stats = "results/figs/ovary_silenced_enrichment/ovary_silenced_enrichment.stats.tsv",
    script:
        "../fig-scripts/silenced-in-ovary.R"

rule collect_computed_statistics:
    input:
        rules.larval_scrna_vs_bulk.output.stats,
        rules.larval_normalized_te_umis.output.stats,
        rules.larval_raw_te_umis.output.stats,
        rules.larval_y_gene_enr.output.stats,
        rules.tidal_larval_tep_xa_boxplot.output.stats,
        rules.larracuente_y_ins_barchart.output.stats,
        rules.larracuente_y_intronic_ins.output.stats,
        rules.w1118_y_linked_copies.output.stats,
        rules.w1118_pct_y_linked_rna_vs_wgs.output.stats,
        rules.w1118_full_length_tx.output.stats,
        rules.comparison_with_mahadevaraju.output.stats,
        rules.all_dataset_tep_scores.output.stats,
        rules.flam_tep_enrichment.output.stats,
        rules.ovary_silenced_enrichment.output.stats
    output:
        stats="results/figs/collect_computed_statistics/collect_computed_statistics.stats.tsv"
    script:
        "../fig-scripts/collect_computed_statistics.R"

rule collect_diffs:
    input:
        diffs = scrna("results/scanpy/{group}/diffs.csv.gz"),
        rename = rules.celltype_rename_table.output.tsv,
    output:
        tsv = "results/finalized/{group}.diffs.tsv.gz"
    resources:
        time=40,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-diffs.R"

rule ica_reproducibility:
    input:
        FINAL_DATASETS_FILES,
    output:
        png1 = "results/figs/ica_reproducibility/ica_reproducibility.pval.png",
        ggp1 = "results/figs/ica_reproducibility/ica_reproducibility.pval.ggp.rds",
        png2 = "results/figs/ica_reproducibility/ica_reproducibility.corr.png",
        ggp2 = "results/figs/ica_reproducibility/ica_reproducibility.corr.ggp.rds",
        dat = "results/figs/ica_reproducibility/ica_reproducibility.dat.tsv",
    script:
        "../fig-scripts/ica-repro.R"

rule spermatogonia_marker_expr:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/spermatogonia_marker_expr/spermatogonia_marker_expr.png",
        ggp = "results/figs/spermatogonia_marker_expr/spermatogonia_marker_expr.ggp.rds",
        dat = "results/figs/spermatogonia_marker_expr/spermatogonia_marker_expr.dat.tsv"
    script:
        "../fig-scripts/spermatogonia-marker-expr.R"

rule te_expression_by_cluster:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/te_expression_by_cluster/te_expression_by_cluster.png",
        ggp = "results/figs/te_expression_by_cluster/te_expression_by_cluster.ggp.rds",
        dat = "results/figs/te_expression_by_cluster/te_expression_by_cluster.dat.tsv"
    script:
        "../fig-scripts/te-expression-by-cluster.R"
