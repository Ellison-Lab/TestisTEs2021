
# ------------------------------------------------------------------------------
# ICA itself
# ------------------------------------------------------------------------------
rule grid_ica_reps:
    input:
        rules.standardize.output
    output:
        source = "results/gep-grid-search/{group}/ica/{k}/{ovr}/{rep}/source.csv.gz",
        mixing = "results/gep-grid-search/{group}/ica/{k}/{ovr}/{rep}/mixing.csv.gz",
    params:
        random_seed = lambda wc: GRID_INDIV_RANDOM_SEEDS_ICA[int(wc.ovr)-1][int(wc.rep)-1],
        comps = lambda wc: int(wc.k)
    resources:
        time=40,
        mem=12000,
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/ica.py"

rule grid_ica_consensus:
    input:
        source= lambda wc: expand("results/gep-grid-search/{g}/ica/{k}/{ovr}/{rep}/source.csv.gz",g=wc.group, k=wc.k,ovr = wc.ovr, rep=range(1,max(REPS)+1,1)),
        mixing= lambda wc: expand("results/gep-grid-search/{g}/ica/{k}/{ovr}/{rep}/mixing.csv.gz",g=wc.group, k=wc.k, ovr = wc.ovr, rep=range(1,max(REPS)+1,1)),
    output:
        usage="results/gep-grid-search/{group}/ica/{k}/{ovr}/consensus-usage.csv.gz",
        ica="results/gep-grid-search/{group}/ica/{k}/{ovr}/consensus-ica.csv.gz",
        dists = "results/gep-grid-search/{group}/ica/{k}/{ovr}/consensus-dists.csv.gz",
        silhouette = "results/gep-grid-search/{group}/ica/{k}/{ovr}/consensus-silhouette.csv.gz",
    params:
        knn = OUTLIER_FILT_KNN_ICA,
        max_dist = OUTLIER_MAX_DIST_ICA,
        k = lambda wc: wc.k
    resources:
        time=240,
        mem=128000,
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/grid-ica-consensus.R"

rule grid_combine_ica_metrics:
    input:
        lambda wc: expand("results/gep-grid-search/{g}/ica/{k}/{ovr}/consensus-{im}.csv.gz",g=wc.group, k=COMPONENTS,ovr=OVERALL_REPS, im=wc.ica_metric)
    output:
        csv = "results/gep-grid-search/{group}/ica/{ica_metric}.csv"
    resources:
        time=60,
        mem=12000,
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/gather_consensus_ica_metrics.R"

# ------------------------------------------------------------------------------
# fdr calcs and cutoffs
# ------------------------------------------------------------------------------

rule grid_fdr_calc:
    input:
        rules.grid_ica_consensus.output.ica,
    output:
        "results/gep-grid-search/{group}/ica/{k}/{ovr}/consensus-ica.qvalues.csv.gz"
    params:
        ICAver = ICA_VERSION
    resources:
        time=120,
        mem=36000,
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/qval_calc.R"

# ------------------------------------------------------------------------------
# enrichment
# ------------------------------------------------------------------------------

rule grid_run_topgo:
    input:
        rules.grid_fdr_calc.output
    output:
        "results/gep-grid-search/{group}/enr/{k}/{ovr}/{fdr}/enr.csv"
    params:
        qval = lambda wc: wc.fdr,
        nodes = TOPGO_NODES,
        ont = ONT
    resources:
        time=240,
        mem=36000,
    conda:
        "../envs/topgo.yaml"
    script:
        "../scripts/go_enr.R"

rule grid_combine_reps_enr:
    input:
        lambda wc: expand("results/gep-grid-search/{g}/enr/{k}/{ovr}/{fdr}/enr.csv",g=wc.group, k=COMPONENTS,ovr=OVERALL_REPS,fdr=QVALS)
    output:
        csv="results/gep-grid-search/{group}/enr.csv"
    resources:
        time=60,
        mem=12000,
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/gather_enr.R"

rule grid_combine_enr_metrics:
    input:
      rules.grid_combine_reps_enr.output
    output:
      "results/gep-grid-search/{group}/enr/enr-metrics.csv"
    params:
        ont = ONT
    conda:
        "../envs/topgo.yaml"
    script:
      "../scripts/gather_supp_enr_metrics.R"
