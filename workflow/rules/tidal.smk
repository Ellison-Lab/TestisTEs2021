rule get_tidal:
    output:
        "results/tidal/download/tidalfly.zip"
    params:
        tidal = config.get("TIDAL")
    shell:
        """
        wget -O {output} {params.tidal}
        """

rule unzip_tidal:
    input:
        rules.get_tidal.output
    output:
        lab = "results/tidal/data/LabStrain_flies.zip",
        dgrp = "results/tidal/data/DGRP_flies.zip",
        dgn = "results/tidal/data/DGN_flies.zip",
        pool = "results/tidal/data/Pool_Flies.zip",
        cells = "results/tidal/data/CellLines.zip",
    conda:
        "../envs/unzip.yaml"
    shell:
        """
        unzip -d results/tidal/data {input}
        """

rule unzip_dgrp:
    input:
        rules.unzip_tidal.output.dgrp
    output:
        directory("results/tidal/data/DGRP")
    conda:
        "../envs/unzip.yaml"
    shell:
        """
        unzip -d {output} {input}
        """

rule get_tidal_anno:
    output:
        "results/tidal/download/annotation.tar.gz"
    params:
        tidal = config.get("TIDAL_ANNO")
    shell:
        """
        wget -O {output} {params.tidal}
        """

rule unzip_tidal_anno:
    input:
        rules.get_tidal_anno.output
    output:
        sizes = "results/tidal/anno/annotation/dm6.chr.len",
        struct = "results/tidal/anno/annotation/fly_virus_structure_repbase.fa",
        mappa = "results/tidal/anno/annotation/gem_mappability_dm6_100mer.mappability",
        reflat = "results/tidal/anno/annotation/refflat_dm6.txt",
        repmask = "results/tidal/anno/annotation/repmasker_dm6_track.txt",
        classes = "results/tidal/anno/annotation/Tidalbase_Dmel_TE_classifications_2015.txt",
        seqs = "results/tidal/anno/annotation/Tidalbase_transposon_sequence.fa",
    params:
        out = "results/tidal/anno"
    shell:
        """
        mkdir -p {params.out} &&
        tar -xzf {input} -C {params.out}
        """

rule tidal_xa_ratio:
    input:
        sizes = rules.unzip_tidal_anno.output.sizes,
        dgrp =  rules. unzip_dgrp.output,
        geps = rules.fdr_cut.output,
        lookup = 'resources/te_id_lookup.curated.tsv.txt',
    output:
        "results/tidal/ratio/{group}/xa-ratio.csv"
    params:
        group = lambda wc: wc.group
    conda:
        "../envs/bioc-general.yaml"
    script:
        "../scripts/tidal-expected-ins.R"
