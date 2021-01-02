CUTOFF_APPROACH_ICA = config.get("CUTOFF_APPROACH_ICA", 1)
#QVAL_ICA = config.get("QVAL_ICA",0.005)
#COMPS_ICA = config.get("COMPS_ICA", 140)
REPS_ICA = config.get("REPS_ICA",10)
INDIV_RANDOM_SEEDS_ICA = random.sample(range(0,100000,1),k=REPS_ICA)

ONTS = config.get("ONTS",["BP","MF","CC"])


def get_optimal_ica(wcg, param='comps'):
    f=open("results/gep-grid-search/{g}/optimal.json".format(g=wcg))
    x = json.load(f)
    return(x[param][0])

rule ica_reps:
    input:
        rules.standardize.output,
        rules.find_optimal_cica_params.output
    params:
        random_seed = lambda wc: INDIV_RANDOM_SEEDS_ICA[int(wc.ica_rep)-1],
        comps = lambda wc: get_optimal_ica(wc.group)
    output:
        source="results/gep/{group}/optimal/reps/{ica_rep}/source.csv.gz",
        mixing="results/gep/{group}/optimal/reps/{ica_rep}/mixing.csv.gz",
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
        source= lambda wc: expand("results/gep/{g}/optimal/reps/{rep}/source.csv.gz",g=wc.group, rep=range(1,REPS_ICA+1,1)),
        mixing= lambda wc: expand("results/gep/{g}/optimal/reps/{rep}/mixing.csv.gz",g=wc.group, rep=range(1,REPS_ICA+1,1)),
    output:
        usage="results/gep/{group}/optimal/consensus-usage.csv.gz",
        ica="results/gep/{group}/optimal/consensus-ica.csv.gz",
        dists = "results/gep/{group}/optimal/consensus-dists.csv.gz",
        silhouette = "results/gep/{group}/optimal/consensus-silhouette.csv.gz",
    params:
        knn = OUTLIER_FILT_KNN_ICA,
        max_dist = OUTLIER_MAX_DIST_ICA,
        k = lambda wc: get_optimal_ica(wc.group)
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
        "results/gep/{group}/optimal/consensus-ica-qvalues.csv.gz"
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
        "results/gep/{group}/optimal/consensus-ica-modules.json"
    params:
        q = lambda wc: get_optimal_ica(wc.group,'qval')
    conda:
        "../envs/gep.yaml"
    script:
        "../scripts/qval_cut.R"

rule enrich_ica_modules:
    input:
        rules.fdr_calc.output
    output:
        "results/gep/{group}/optimal/consensus-ica-enrichment-{ont}.csv.gz"
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
