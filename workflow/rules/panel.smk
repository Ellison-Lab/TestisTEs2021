rule figure1:
    input:
        rules.intro_larval_umap.output,
        rules.larval_marker_expression.output,
        rules.larval_scrna_basic_qc_stats.output
    output:
        "results/panels/figure1.pdf",
        "results/panels/figure1.rds"
    script:
        "../fig-scripts/figure1.R"

rule figure2:
    input:
        rules.larval_te_heatmap.output,
        rules.te_expression_by_cluster.output
    output:
        "results/panels/figure2.pdf",
        "results/panels/figure2.expr-in-half.pdf",
        "results/panels/figure2.rds"
    script:
        "../fig-scripts/figure2.R"

rule figure3:
    input:
        rules.larval_te_heatmap.output,
        rules.larval_tep_pie.output,
        rules.larval_tep_te_barchart.output,
        rules.larval_all_gep_te_bachart.output,
        rules.larval_tep_usage_umap.output,
    output:
        "results/panels/figure3.pdf",
        "results/panels/figure3.rds"
    script:
        "../fig-scripts/figure3.R"

rule figure4:
    input:
        rules.larval_later_sperm_marker_umis.output,
        rules.larval_y_gene_dotplot.output,
        rules.larval_fish_candidate_umis.output,
    output:
        "results/panels/figure4.pdf",
        "results/panels/figure4.rds"
    script:
        "../fig-scripts/figure4.R"

rule figure5:
    input:
        rules.tidal_larval_tep_xa_boxplot.output,
        rules.larracuente_y_ins_barchart.output,
        rules.w1118_pct_y_linked_rna_vs_wgs.output,
        rules.w1118_y_linked_copies.output,
        rules.larval_pirna_genes.output,
    output:
        "results/panels/figure5.pdf",
        "results/panels/figure5.rds"
    script:
        "../fig-scripts/figure5.R"

rule supp1:
    input:
        rules.comparison_with_mahadevaraju.output,
        rules.larval_scrna_vs_bulk.output,
        rules.spermatogonia_marker_expr.output,
        rules.larval_scrna_basic_qc_stats.output
    output:
        "results/panels/supp1.pdf",
        "results/panels/supp1.rds"
    script:
        "../fig-scripts/supp1.R"

rule supp2:
    input:
        rules.w1118_full_length_tx.output,
        rules.chimerics.output,
        rules.larval_raw_te_umis.output,
        rules.larval_normalized_te_umis.output,
    output:
        "results/panels/supp2.pdf",
        "results/panels/supp2.rds"
    script:
        "../fig-scripts/supp2.R"

rule supp3:
    input:
        rules.larval_cica_grid_heatmap.output,
        rules.ica_reproducibility.output,
        rules.larval_gep_size_histogram.output,
        rules.all_dataset_tep_scores.output,
        rules.larval_ica_optimized_enr.output
    output:
        "results/panels/supp3.pdf",
        "results/panels/supp3.rds"
    script:
        "../fig-scripts/supp3.R"

rule supp4:
    input:
        rules.larval_gep_corr_heatmap.output,
        rules.larval_tep_pie.output,
        rules.larval_tep_te_barchart.output
    output:
        "results/panels/supp4.pdf",
        "results/panels/supp4.rds"
    script:
        "../fig-scripts/supp4.R"

rule supp4_5:
    input:
        rules.flam_tep_enrichment.output,
        rules.ovary_silenced_enrichment.output
    output:
        "results/panels/supp4_5.pdf",
        "results/panels/supp4_5.rds"
    script:
        "../fig-scripts/supp4-5.R"

rule supp5:
    input:
        rules.larval_later_sperm_marker_umis.output,
        rules.larval_y_gene_enr.output,
        rules.larval_y_gene_dotplot.output,
        rules.w1118_y_linked_copies.output,
        rules.larval_tep_go_enrichment.output,
    output:
        "results/panels/supp5.pdf",
        "results/panels/supp5.rds"
    script:
        "../fig-scripts/supp5.R"

rule supp_figs:
    input:
        rules.figure1.output,
        rules.figure2.output,
        rules.figure3.output,
        rules.figure4.output,
        rules.figure5.output,
        rules.supp1.output,
        rules.supp2.output,
        rules.supp3.output,
        rules.supp4.output,
        rules.supp4_5.output,
        rules.supp5.output,
        rules.larracuente_y_intronic_ins.output,
        rules.arriba_total_rna_fusions.output,
        #rules.collect_computed_statistics.output,
    output:
        "results/collected/supp.pdf"
    script:
        "../markdown/figures.Rmd"

rule main_figs:
    input:
        rules.figure1.output,
        rules.figure2.output,
        rules.figure3.output,
        rules.figure4.output,
        rules.figure5.output,
        rules.supp1.output,
        rules.supp2.output,
        rules.supp3.output,
        rules.supp4.output,
        rules.supp4_5.output,
        rules.supp5.output,
        rules.larracuente_y_intronic_ins.output,
        rules.arriba_total_rna_fusions.output,
        #rules.collect_computed_statistics.output,
    output:
        f1 = "results/collected/figure1.pdf",
        f2 = "results/collected/figure2.pdf",
        f3 = "results/collected/figure3.pdf",
        f4 = "results/collected/figure4.pdf",
        f5 = "results/collected/figure5.pdf"
    script:
        "../fig-scripts/figures-pg.R"

rule get_source_data:
    input:
        rules.main_figs.output,
        rules.supp_figs.output
    output:
        "results/collected/SourceData.xlsx"
    script:
        "../scripts/collect_source_data.R"
