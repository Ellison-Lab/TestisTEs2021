CUTOFF_APPROACH_ICA = config.get("CUTOFF_APPROACH_ICA", 1)
QVAL_ICA = config.get("QVAL_ICA",0.005)
COMPS_ICA = config.get("COMPS_ICA", 140)
REPS_ICA = config.get("REPS_ICA",10)
INDIV_RANDOM_SEEDS_ICA = random.sample(range(0,100000,1),k=REPS_ICA)
OUTLIER_FILT_KNN_ICA = config.get("OUTLIER_FILT_KNN_ICA",5)
OUTLIER_MAX_DIST_ICA = config.get("OUTLIER_MAX_DIST_ICA",250)
MAX_SCALE_ICA = config.get("MAX_SCALE_ICA")

ONTS = config.get("ONTS",["BP","MF","CC"])

rule standardize:
    input:
        expr = "results/scanpy/{group}/lognorm-expression.csv.gz",
        genes = "results/scanpy/{group}/hivar-and-tes.csv",
    output:
        "results/gep/{group}/standardized.csv.gz"
    params:
        maxval = MAX_SCALE_ICA
    resources:
        time=40,
        mem=12000,
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/standardize.py"

rule ica_reps:
    input:
        rules.standardize.output,
    params:
        random_seed = lambda wc: INDIV_RANDOM_SEEDS_ICA[int(wc.ica_rep)-1],
        comps = lambda wc: wc.k
    output:
        source="results/gep/{group}/{k}/reps/{ica_rep}/source.csv.gz",
        mixing="results/gep/{group}/{k}/reps/{ica_rep}/mixing.csv.gz",
        #components="out/ica/{k}/reps/{ica_rep}/components.csv.gz",
    resources:
        time=30,
        mem=24000,
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/ica.py"

rule ica_consensus:
    input:
        source= lambda wc: expand("results/gep/{g}/{k}/reps/{rep}/source.csv.gz",g=wc.group, k=wc.k,rep=range(1,REPS_ICA+1,1)),
        mixing= lambda wc: expand("results/gep/{g}/{k}/reps/{rep}/mixing.csv.gz",g=wc.group, k=wc.k, rep=range(1,REPS_ICA+1,1)),
    output:
        usage="results/gep/{group}/{k}/consensus-usage.csv.gz",
        ica="results/gep/{group}/{k}/consensus-ica.csv.gz",
        dists = "results/gep/{group}/{k}/consensus-dists.csv.gz",
        silhouette = "results/gep/{group}/{k}/consensus-silhouette.csv.gz",
    params:
        knn = OUTLIER_FILT_KNN_ICA,
        max_dist = OUTLIER_MAX_DIST_ICA,
        k = COMPS_ICA
    resources:
        time=120,
        mem=24000,
        cpus=2
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/ica-consensus.R"

rule fdr_calc:
    input:
        rules.ica_consensus.output.ica
    output:
        "results/gep/{group}/{k}/consensus-ica-qvalues.csv.gz"
    params:
        ICAver = CUTOFF_APPROACH_ICA
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/qval_calc.R"

rule fdr_cut:
    input:
        rules.fdr_calc.output
    output:
        "results/gep/{group}/{k}/consensus-ica-modules.json"
    params:
        q = QVAL_ICA
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/qval_cut.R"

rule enrich_ica_modules:
    input:
        rules.fdr_calc.output
    output:
        "results/gep/{group}/{k}/consensus-ica-enrichment-{ont}.csv.gz"
    params:
        qval = QVAL_ICA,
        ont = lambda wc: wc.ont,
        nodes = 100,
    resources:
        time=60,
        mem=12000,
        cpus=2
    conda:
        "../envs/topgo.yaml"
    script:
        "../scripts/enr.R"
