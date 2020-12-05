rule remove_hashtags_from_te_fa:
    input:
        a = determine_resource(config.get("TRANSPOSON_FASTA"))
    output:
        "results/custom-genome/consensus_te.fasta"
    conda:
        "../envs/bioc-general.yaml"
    group:
        "custom-genome"
    log:
        "results/logs/remove_hashtags_from_te_fa/remove_hashtags_from_te_fa.log"
    script:
        "../scripts/remove_hashtags_from_te_fa.R"

rule idx_consensus_te_fa:
    input:
        rules.remove_hashtags_from_te_fa.output
    output:
        "results/custom-genome/consensus_te.fasta.fai"
    conda:
        "../envs/samtools.yaml"
    group:
        "custom-genome"
    log:
        "results/logs/idx_consensus_te_fa/idx_consensus_te_fa.log"
    shell:
        """
        samtools faidx {input} 2> {log}
        """

rule consensus_te_fa_2_gtf:
    input:
        rules.idx_consensus_te_fa.output
    conda:
        "../envs/bioc-general.yaml"
    output:
        temp("results/custom-genome/consensus_te_fa.gtf")
    group:
        "custom-genome"
    log:
        "results/logs/consensus_te_fa_2_gtf/consensus_te_fa_2_gtf.log"
    script:
        "../scripts/ashburner_consensus_te_fai_2_gtf.R"

rule combine_fas:
    input:
        tes = rules.remove_hashtags_from_te_fa.output,
        g   = determine_resource(config.get("GENOME_FASTA"))
    output:
        "results/custom-genome/combined.fasta",
    group:
        "custom-genome"
    log:
        "results/logs/combine_fas/combine_fas.log"
    shell:
        """
        cat {input.tes} {input.g} > {output} 2> {log}
        """

rule combine_gtfs:
    input:
        tes = rules.consensus_te_fa_2_gtf.output,
        genes = determine_resource(config.get("GENE_GTF"))
    output:
        temp("results/custom-genome/combined.gtf")
    group:
        "custom-genome"
    log:
        "results/logs/combine_gtfs/combine_gtfs.log"
    shell:
        """
        cat <(zcat {input.genes}) {input.tes} > {output} 2> {log}
        """

rule fix_strand_gtf:
    input:
        rules.combine_gtfs.output
    output:
        "results/custom-genome/combined.fixed.gtf"
    params:
        disallow = config.get("DISALLOW_GENE_TYPES")
    conda:
        "../envs/bioc-general.yaml"
    group:
        "custom-genome"
    log:
        "results/logs/fix_strand_gtf/fix_strand_gtf.log"
    script:
        "../scripts/fix_strand_combined_gtfs.R"
