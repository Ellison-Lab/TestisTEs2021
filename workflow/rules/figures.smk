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
        dat = "results/figs/larval_normalized_te_umis/larval_normalized_te_umis.dat.tsv"
    script:
        "../fig-scripts/larval_normalized_te_umis.R"



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


rule larval_tep_te_bachart:
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

rule larval_gep_corr_heatmap:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png1 = "results/figs/larval_gep_corr_heatmap/larval_gep_corr.usage.png",
        png2 = "results/figs/larval_gep_corr_heatmap/larval_gep_corr.membership.png",
    script:
        "../fig-scripts/larval_gep_corr_heatmap.R"

rule larval_y_gene_dotplot:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.png",
        ggp = "results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.ggp.rds",
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
        png = "results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.png",
        ggp = "results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.ggp.rds",
        dat = "results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.dat.tsv"
    script:
        "../fig-scripts/tidal_larval_tep_xa_boxplot.R"

rule larracuente_y_ins_barchart:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.png",
        ggp = "results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.ggp.rds",
        dat = "results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.dat.tsv"
    script:
        "../fig-scripts/larracuente_y_ins_barchart.R"


rule w1118_y_linked_copies:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.png",
        ggp = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.ggp.rds",
        dat = "results/figs/w1118_y_linked_copies/w1118_y_linked_copies.dat.tsv"
    script:
        "../fig-scripts/w1118_y_linked_copies.R"

rule w1118_pct_y_linked_rna_vs_wgs:
    input:
        FINAL_DATASETS_FILES,
        rules.celltype_rename_table.output
    output:
        png = "results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.png",
        ggp = "results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.ggp.rds",
        dat = "results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.dat.tsv"
    script:
        "../fig-scripts/w1118_pct_y_linked_rna_vs_wgs.R"

# ------------------------------------------------------------------------------

# rule larval_pirna_expression:
#     input:
#         FINAL_DATASETS_FILES,
#         rules.celltype_rename_table.output
#     output:
#         png = "results/figs/larval_pirna_expression/larval_pirna_expression.png",
#         ggp = "results/figs/larval_pirna_expression/larval_pirna_expression.ggp.rds",
#         dat = "results/figs/larval_pirna_expression/larval_pirna_expression.dat.tsv"
#     script:
#         "../fig-scripts/larval_pirna_expression.R"
#
# rule larval_chromatin_mod_umap:
#     input:
#         FINAL_DATASETS_FILES,
#         rules.celltype_rename_table.output
#     output:
#         png = "results/figs/larval_chromatin_mod_umap/larval_chromatin_mod_umap.png",
#         ggp = "results/figs/larval_chromatin_mod_umap/larval_chromatin_mod_umap.ggp.rds",
#         dat = "results/figs/larval_chromatin_mod_umap/larval_chromatin_mod_umap.dat.tsv"
#     script:
#         "../fig-scripts/larval_chromatin_mod_umap.R"
#
#
#
# rule larval_chrom_modifiers_umis:
#     input:
#         FINAL_DATASETS_FILES,
#         rules.celltype_rename_table.output
#     output:
#         png = "results/figs/larval_chrom_modifiers_umis/larval_chrom_modifiers_umis.png",
#         ggp = "results/figs/larval_chrom_modifiers_umis/larval_chrom_modifiers_umis.ggp.rds",
#         dat = "results/figs/larval_chrom_modifiers_umis/larval_chrom_modifiers_umis.dat.tsv"
#     script:
#         "../fig-scripts/larval_chrom_modifiers_umis.R"
#
# rule comparison_with_mahadevaraju:
#     input:
#         FINAL_DATASETS_FILES,
#         rules.celltype_rename_table.output
#     output:
#         png = "results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.png",
#         ggp = "results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.ggp.rds",
#         dat = "results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.dat.tsv"
#     script:
#         "../fig-scripts/comparison_with_mahadevaraju.R"
