"""
Link to doc with dl link for pacbio barcodes
https://www.pacb.com/wp-content/uploads/Sequel_RSII_96_barcodes_v1.fasta_.zip
https://www.pacb.com/wp-content/uploads/Procedure-Checklist-Preparing-SMRTbell-Libraries-using-PacBio-Barcoded-Universal-Primers-for-Multiplexing-Amplicons.pdf
"""

rule flnc_bam_from_fq:
    input:
        rules.fastq_dump.output
    output:
        "results/isoseq-1/flnc/{sample}.bam"
    conda:
        "../envs/picard.yaml"
    resources:
        mem=16000,
    params:
        reads = lambda wc: "results/fastq-dump/{s}/{r}.fastq".format(s=wc.sample, r= basename(pep.get_sample(wc.sample).sra_accession[0])),
    shell:
        """
        picard -Xmx{resources.mem}m FastqToSam --FASTQ {params.reads} -O {output} --SAMPLE_NAME {wildcards.sample}
        """

rule isoseq_cluster:
    """
    https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md
    """
    input:
        reads = rules.flnc_bam_from_fq.output
    output:
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.bam",
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.hq.fasta.gz",
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.lq.fasta.gz",
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.hq.bam",
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.lq.bam",
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.transcriptset.xml",
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.cluster",
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.cluster_report.csv",
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.bam.pbi",
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.hq.bam.pbi",
        "results/isoseq-1/isoseq_cluster/{sample}/{sample}.lq.bam.pbi",
    conda:
        "../envs/isoseq3.yaml"
    log:
        "results/logs/isoseq_cluster/{sample}.log"
    threads:
        24
    resources:
        time=240,
        mem=128000,
        cpus=24
    shell:
        """
        isoseq3 cluster -j {threads} --use-qvs --verbose {input} {output[0]} 2> {log}
        """

rule pbmm2_isoseq:
    input:
        txs = rules.isoseq_cluster.output[5],
        ref = rules.combine_fas.output,
    output:
        "results/isoseq-1/pbmm2/{sample}/{sample}.bam",
    conda:
        "../envs/minimap2.yaml"
    log:
        "results/logs/pbmm2_isoseq/{sample}.log"
    threads:
        32
    resources:
        time=60,
        mem=32000,
        cpus=32
    shell:
        """
        pbmm2 align -j {threads} --sort --log-file {log} --preset ISOSEQ {input.ref} {input.txs} {output[0]}
        """




rule isoseq_collapse:
    input:
        reads = rules.pbmm2_isoseq.output[0]
    output:
        "results/isoseq-1/isoseq_collapse/{sample}.gff"
    conda:
        "../envs/isoseq3.yaml"
    log:
        "results/logs/isoseq_collapse/{sample}.log"
    threads:
        12
    resources:
        time=60,
        mem=16000,
        cpus=12
    shell:
        """
        isoseq3 collapse -j {threads} {input} {output}
        """

rule mm2_isoseq_fusion:
    input:
        reads = rules.isoseq_cluster.output[1],
        ref = rules.combine_fas.output,
    output:
        "results/isoseq-1/mm2_isoseq_fusion/{sample}/{sample}.sam"
    conda:
        "../envs/minimap2.yaml"
    log:
        "results/logs/mm2_isoseq_fusion/{sample}.log"
    threads:
        30
    resources:
        time=30,
        mem=32000,
        cpus=32
    shell:
        """
        minimap2 -ax splice:hq -t {threads} --secondary=no {input.ref} {input.reads} | \
            samtools sort -O SAM > {output}
        """
rule gunzip_pre_cupcake_fusion:
    input:
        reads = rules.isoseq_cluster.output[1],
    output:
        "results/isoseq-1/cupcake-fusion/{sample}/{sample}.hq.fasta",
    shell:
        "gunzip -c {input} > {output}"

rule cupcake_fusion:
    input:
        sam = rules.mm2_isoseq_fusion.output,
        reads = rules.gunzip_pre_cupcake_fusion.output,
        report = rules.isoseq_cluster.output[7]
    output:
        "results/isoseq-1/cupcake-fusion/{sample}/{sample}.fusion.gff",
    params:
        pfx = lambda wc: "results/isoseq-1/cupcake-fusion/{s}/{s}.fusion".format(s=wc.sample),
    conda:
        "../envs/cupcake.yaml"
    resources:
        time=30,
        mem=16000,
        cpus=1
    log:
        "results/logs/cupcake_fusion/{sample}.log"
    shell:
        """
        fusion_finder.py --input {input.reads} -s {input.reads} -o {params.pfx} \
            --cluster_report_csv {input.report}
            2> {log}
        """


rule mm2_isoseq:
    input:
        reads = rules.fastq_dump.output,
        ref = rules.combine_fas.output,
    output:
        "results/isoseq-1/mm2/{sample}.sam"
    conda:
        "../envs/minimap2.yaml"
    log:
        "results/logs/mm2_isoseq/{sample}.log"
    params:
        reads = lambda wc: "results/fastq-dump/{s}/{r}.fastq".format(s=wc.sample, r= basename(pep.get_sample(wc.sample).sra_accession[0])),
    threads:
        30
    group:
        "isoseq-1-map"
    resources:
        time=30,
        mem=32000,
        cpus=32
    shell:
        """
        minimap2 -ax splice:hq -t {threads} --secondary=no {input.ref} {params.reads} | \
            samtools sort -O SAM > {output}
        """

rule mm2_isoseq_2_bam:
    input:
        rules.mm2_isoseq.output
    output:
        "results/isoseq-1/mm2/{sample}.bam"
    conda:
        "../envs/minimap2.yaml"
    group:
        "isoseq-1-map"
    log:
        "results/logs/mm2_isoseq_bai/{sample}.log"
    shell:
        """
        samtools view -O BAM {input} > {output}
        """

rule mm2_isoseq_bai:
    input:
        rules.mm2_isoseq_2_bam.output
    output:
        "results/isoseq-1/mm2/{sample}.bam.bai"
    conda:
        "../envs/minimap2.yaml"
    group:
        "isoseq-1-map"
    log:
        "results/logs/mm2_isoseq_bai/{sample}.log"
    shell:
        """
        samtools index {input}
        """


rule cupcake_collapse:
    input:
        sam = rules.mm2_isoseq.output,
        reads = rules.fastq_dump.output,
    output:
        "results/isoseq-1/cupcake/{sample}/{sample}.collapsed.gff",
        "results/isoseq-1/cupcake/{sample}/{sample}.collapsed.rep.fq",
        "results/isoseq-1/cupcake/{sample}/{sample}.collapsed.group.txt,"
    params:
        reads = lambda wc: "results/fastq-dump/{s}/{r}.fastq".format(s=wc.sample, r= basename(pep.get_sample(wc.sample).sra_accession[0])),
    conda:
        "../envs/cupcake.yaml"
    resources:
        time=30,
        mem=32000,
        cpus=1
    log:
        "results/logs/cupcake_collapse/{sample}.log"
    shell:
        """
        collapse_isoforms_by_sam.py --input {params.reads} --fq \
            -s {input.sam} --dun-merge-5-shorter \
            -o results/isoseq-1/cupcake/{wildcards.sample}/{wildcards.sample} \
            2> {log}
        """
