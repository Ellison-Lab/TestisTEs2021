rule collect_var:
    input:
        scrna("results/scanpy/{group}/var.csv/")
    output:
        directory("results/finalized/{group}/var")
    resources:
        time=40,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_obs:
    input:
        clusters = scrna("results/scanpy/{group}/clusters-to-celltypes.csv"),
        obsm = scrna("results/scanpy/{group}/obsm.csv"),
        obs = scrna("results/scanpy/{group}/obs.csv"),
    output:
        directory("results/finalized/{group}/obs")
    resources:
        time=40,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-obs.R"

rule collect_scrna_expr:
    input:
        var= rules.collect_var.output,
        expr = scrna("results/scanpy/{group}/lognorm-expression.csv.gz"),
        obs = rules.collect_obs.output
    output:
        directory("results/finalized/{group}/expr")
    resources:
        time=40,
        mem=128000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-scrna-expr.R"

rule collect_scrna_scaled:
    input:
        var= rules.collect_var.output,
        expr = scrna("results/scanpy/{group}/scaled-expression.csv.gz"),
        obs = rules.collect_obs.output
    output:
        directory("results/finalized/{group}/scaled")
    resources:
        time=40,
        mem=128000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-scrna-expr.R"

rule collect_umis:
    input:
        scrna("results/scanpy/{group}/umis.csv.gz")
    output:
        directory("results/finalized/{group}/umis")
    resources:
        time=40,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_grid_enr_metrics:
    input:
        gep("results/grid-search-{group}-go-metrics.csv")
    output:
        directory("results/finalized/{group}/grid_enr")
    resources:
        time=10,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

# rule collect_grid_cica_silhouette:
#     input:
#         "results/gep-grid-search/{group}/silhouette.csv"
#     output:
#         directory("results/finalized/{group}/grid_silhouette")
#     resources:
#         time=10,
#         mem=24000,
#     conda:
#         "../envs/r_arrow.yaml"
#     script:
#         "../scripts/collect-var.R"

rule collect_optimal_enr:
    input:
        lambda wc: gep(expand("results/gep/{g}/optimal/consensus-ica-enrichment-{o}.csv.gz",g = wc.group, o=config.get('ONTS')))
    output:
        directory("results/finalized/{group}/optimal_gep_enr/")
    resources:
        time=10,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-gep-enr.R"

rule collect_optimal_usage:
    input:
        lambda wc: gep(expand("results/gep/{g}/optimal/consensus-usage.csv.gz",g = wc.group, o=config.get('ONTS')))
    output:
        directory("results/finalized/{group}/optimal_gep_usage/")
    resources:
        time=10,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_optimal_gep_membership:
    input:
        lambda wc: gep(expand("results/gep/{g}/optimal/consensus-ica-qvalues.csv.gz",g = wc.group, o=config.get('ONTS')))
    output:
        directory("results/finalized/{group}/optimal_gep_membership/")
    resources:
        time=10,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_optimal_gep_params:
    input:
        gep("results/grid-search-{group}-optimal.json")
    output:
        "results/finalized/optimal-gep-params/{group}.json"
    shell:
        "cp {input} {output}"

rule collect_hetchrom_assembly_ins:
    input:
        hetchrom_y("results/hetchrom-ins/insertions.csv")
    output:
        directory("results/finalized/hetchrom_assembly_insertions/")
    resources:
        time=10,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_tidal:
    input:
        xa_ratio("results/tidal/ratio/{group}/xa-ratio.csv")
    output:
        directory("results/finalized/{group}/xa_ratio/")
    resources:
        time=10,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_te_copies:
    input:
        transposon_variants('results/copies/{wgs_group}.csv'),
    output:
        directory("results/finalized/wgs/{wgs_group}/copies/")
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_male_te_wgs_snp_depth:
    input:
        transposon_variants("results/snps/depth-at-snps.csv.gz"),
    output:
        directory("results/finalized/wgs/w1118_male/snp_depths/")
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_te_expr_snp_depth:
    input:
        te_var_expr("results/depths/w1118_testes-{total_rna_rep}-{total_rna_res_type}-at-male-snps.csv.gz")
    output:
        directory("results/finalized/w1118-testes-total-rna/{total_rna_rep}-{total_rna_res_type}-at-male-snps/")
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule copy_total_rna_te_bws:
    input:
        te_var_expr(expand("results/bigwigs/tes/{s}.{sub}.tes.strand-{t}.rpkm.bw",s=['w1118_testes'],sub=['rep1','rep2','rep3','rep4'],t=['forward','reverse']))
    output:
        expand("results/finalized/bigwigs/total-rna/{s}.{sub}.tes.strand-{t}.rpkm.bw",s=['w1118_testes'],sub=['rep1','rep2','rep3','rep4'],t=['forward','reverse'])
    shell:
        """
        cp {input} results/finalized/bigwigs/total-rna/
        """

rule copy_polya_rna_bws:
    input:
        larval_polya(expand("results/bigwigs/{s}.strand-{t}.rpkm.bw",s=['larval_testes_cleaned_papain_01','larval_testes_papain_02','larval_testes_papain_03','larval_testes_papain_04'],t=['forward','reverse']))
    output:
        expand("results/finalized/bigwigs/polya-rna/{s}.strand-{t}.rpkm.bw",s=['larval_testes_cleaned_papain_01','larval_testes_papain_02','larval_testes_papain_03','larval_testes_papain_04'],t=['forward','reverse'])
    shell:
        """
        cp {input} results/finalized/bigwigs/polya-rna/
        """


rule collect_larval_polya_expr:
    input:
        larval_polya("results/star/{larval_polya_sample}/")
    output:
        "results/finalized/larval-polya/{larval_polya_sample}.tsv"
    shell:
        """
        echo 'gene_id\ttotal\tfirst_strand\tsecond_strand' > {output}
        grep -v 'N_' {input}/ReadsPerGene.out.tab >> {output}
        """

rule collect_wgs_pileups:
    input:
        transposon_variants("results/pileups/w1118_male.pileups.csv.gz"),
    output:
        directory("results/finalized/wgs/w1118_male/pileups/")
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_male_snp_vcf:
    input:
        transposon_variants('results/snps/snps.vcf'),
    output:
        "results/finalized/wgs/w1118_male/snps.vcf"
    shell:
        """
        cp {input} {output}
        """

rule collect_mod_scores:
    input:
        "results/finalized/larval-w1118-testes/optimal_gep_membership/",
        "results/finalized/optimal-gep-params/larval-w1118-testes.json"
    output:
        "results/finalized/x-dataset-comparison/mod_scores.csv.gz",
        "results/finalized/x-dataset-comparison/te_expression.csv.gz"
    script:
        "../scripts/collect_mod_scores.py"

rule collect_pirna_kd_rnaseq:
    input:
        pirna_kd_rnaseq("results/deseq2/pirna_kd_vs_control.res.tsv")
    output:
        "results/finalized/pirna_kd_rnaseq/pirna_kd_vs_control.res.tsv"
    shell:
        """
        cp {input} {output}
        """

rule collect_total_rna_fusions:
    input:
        total_rna_fusions("results/arriba/w1118_testes_{total_rna_fusions_sample}/fusions.tsv")
    output:
        "results/finalized/arriba/w1118_testes_{total_rna_fusions_sample}.fusions.tsv"
    shell:
        """
        cp {input} {output}
        """
