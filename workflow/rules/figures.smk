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
        png = "results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.png",
        ggp = "results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.ggp.rds",
        dat = "results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.dat.tsv",
    script:
        "../fig-scripts/larval_scrna_vs_bulk.R"

rule larval_te_heatmap:


rule larval_normalized_te_umis:

rule larval_raw_te_umis:

rule larval_tep_usage_umap:

rule larval_tep_te_bachart:

rule larval_tep_pie:

rule larval_gep_usage_heatmap:

rule larval_gep_membership_heatmap:

rule larval_cica_grid_heatmap:

rule larval_cica_k_vs_score:

rule larval_gep_size_histogram:

rule larval_gep_corr_by_membership_heatmap:

rule larval_gep_usage_by_membership_heatmap:

rule larval_y_gene_dotplot:

rule larval_gep_term_table:

rule larval_chromatin_mod_umap:

rule larval_chrom_modifiers:

rule larval_fish_candidate_umis:

rule tidal_larval_tep_xa_boxplot:

rule larval_later_sperm_marker_umis:

rule larval_scrna_vs_bulk:

rule comparison_with_mahadevaraju:

rule larracuente_y_ins_barchart:

rule w1118_y_linked_copies:

rule larval_pirna_expression:

rule w1118_pct_y_linked_rna_vs_wgs:
