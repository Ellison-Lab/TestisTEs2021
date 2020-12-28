rule collect_var:
    input:
        "results/scanpy/{group}/anno/var.csv"
    output:
        directory("results/finalized/{group}/var")
    resources:
        time=40,
        mem=12000,
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
        mem=12000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-obs.R"

rule collect_scrna_expr:
    input:
        var= "results/scanpy/{group}/anno/var.csv",
        expr = "results/scanpy/{group}/lognorm-expression.csv.gz"
    output:
        directory("results/finalized/{group}/expr")
    resources:
        time=40,
        mem=12000,
    conda:
        "../envs/r_arrow.yaml"
    script:
        "../scripts/collect-scrna-expr.R"
