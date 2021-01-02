rule collect_var:
    input:
        "results/scanpy/{group}/anno/var.csv"
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
        clusters = "results/scanpy/{group}/clusters-to-celltypes.csv",
        obsm = "results/scanpy/{group}/anno/obsm.csv",
        obs = "results/scanpy/{group}/anno/obs.csv",
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
        expr = "results/scanpy/{group}/lognorm-expression.csv.gz",
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
        expr = "results/scanpy/{group}/scaled-expression.csv.gz",
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

rule collect_grid_enr_metrics:
    input:
        "results/gep-grid-search/{group}/enr-metrics.csv"
    output:
        directory("results/finalized/{group}/grid_enr")
    resources:
        time=10,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_grid_cica_silhouette:
    input:
        "results/gep-grid-search/{group}/silhouette.csv"
    output:
        directory("results/finalized/{group}/grid_silhouette")
    resources:
        time=10,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"

rule collect_optimal_enr:
    input:
        lambda wc: expand("results/gep/{g}/optimal/consensus-ica-enrichment-{o}.csv.gz",g = wc.group, o=config.get('ONTS'))
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
        lambda wc: expand("results/gep/{g}/optimal/consensus-usage.csv.gz",g = wc.group, o=config.get('ONTS'))
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
        lambda wc: expand("results/gep/{g}/optimal/consensus-ica-qvalues.csv.gz",g = wc.group, o=config.get('ONTS'))
    output:
        directory("results/finalized/{group}/optimal_gep_membership/")
    resources:
        time=10,
        mem=24000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-var.R"
