rule chip_trim_se:
    """
    Trim single-end fastqs with fastp.
    """
    input:
        fq = rules.fastq_dump.output,
        fq_input = rules.fastq_dump_input.output,
    output:
        r1 = "results/trimmed/{sample}/{sample}_r1.trimmed.fq.gz",
        html = "results/trimmed/{sample}/{sample}_fastp.html",
        json = "results/trimmed/{sample}/{sample}_fastp.json",
        r1_input = "results/trimmed/{sample}/{sample}_input_r1.trimmed.fq.gz",
        html_input = "results/trimmed/{sample}/{sample}_input_fastp.html",
        json_input = "results/trimmed/{sample}/{sample}_input_fastp.json"
    resources:
        mem=8000,
        cpus=8,
        time=20,
    threads:
        8
    log:
        "results/logs/chip_trim_se/{sample}.log"
    params:
        reads = lambda wc: "results/fastq-dump/{s}/{r}.fastq".format(s=wc.sample, r= basename(pep.get_sample(wc.sample).sra_accession[0])),
        reads_input = lambda wc: "results/fastq-dump/{s}-input/{r}.fastq".format(s=wc.sample, r= basename(pep.get_sample(wc.sample).sra_accession_input[0])),
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp --in1 {params.reads} --out1 {output.r1} \
            -j {output.json} -h {output.html} \
            -w {threads} -L &&

        fastp --in1 {params.reads_input} --out1 {output.r1_input} \
            -j {output.json_input} -h {output.html_input} \
            -w {threads} -L
        """

rule bt2_index_combined:
    input:
        rules.combine_fas.output
    output:
        directory("results/bowtie2-idx/combined")
    resources:
        mem=8000,
        cpus=8,
        time=20,
    threads:
        8
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p {output} &&
        bowtie2-build --threads {threads} {input} {output}/idx
        """

rule aln_chipseq:
    input:
        chip = rules.chip_trim_se.output.r1,
        idx = rules.bt2_index_combined.output
    output:
        "results/chipseq/aln/{sample}/{sample}.bam"
    conda:
        "../envs/bowtie2.yaml"
    resources:
        mem=24000,
        cpus=24,
        time=60,
    threads:
        22
    shell:
        """
        bowtie2 --phred33 -p {threads} \
            --no-discordant --no-unal -k 1 \
            -x {input.idx}/idx {input.chip} | \
            samtools sort -O BAM > {output}
        """

rule aln_chipseq_input:
    input:
        input = rules.chip_trim_se.output.r1_input,
        idx = rules.bt2_index_combined.output
    output:
        "results/chipseq/aln/{sample}/{sample}_input.bam"
    conda:
        "../envs/bowtie2.yaml"
    resources:
        mem=24000,
        cpus=24,
        time=60,
    threads:
        22
    shell:
        """
        bowtie2 --phred33 -p {threads} \
            --no-discordant --no-unal -k 1 \
            -x {input.idx} {input.input} | \
            samtools sort -O BAM > {output}
        """
