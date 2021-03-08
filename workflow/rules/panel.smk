rule figure1:
    input:
        rules.intro_larval_umap.output,
        rules.larval_marker_expression.output
    output:
        "results/panels/figure1.pdf"
    script:
        "../fig-scripts/figure1.R"

rule figure2:
    input:
        rules.larval_te_heatmap.output,
    output:
        "results/panels/figure2.pdf"
    script:
        "../fig-scripts/figure2.R"

rule figure3:
    input:
        rules.larval_te_heatmap.output,
        rules.larval_tep_pie.output,
        rules.larval_tep_te_barchart.output,
    output:
        "results/panels/figure3.pdf"
    script:
        "../fig-scripts/figure3.R"

rule figure4:
    input:
        rules.larval_later_sperm_marker_umis.output,
        rules.larval_y_gene_dotplot.output
    output:
        "results/panels/figure4.pdf"
    script:
        "../fig-scripts/figure4.R"

rule figure5:
    input:
        rules.tidal_larval_tep_xa_boxplot.output,
        rules.larracuente_y_ins_barchart.output,
        rules.w1118_pct_y_linked_rna_vs_wgs.output,
        rules.w1118_y_linked_copies.output,
    output:
        "results/panels/figure5.pdf"
    script:
        "../fig-scripts/figure5.R"
