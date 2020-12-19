rule get_sra:
    params:
        url = lambda wc: pep.get_sample(wc.sample).sra_accession,
        dir = "results/get-sra/"
    output:
        temp(directory("results/get-sra/{sample}"))
    shell:
        """
        wget -P {params.dir}/{wildcards.sample} {params.url}
        """

rule get_sra_input:
    params:
        url = lambda wc: pep.get_sample(wc.sample).sra_accession_input,
        dir = "results/get-sra/"
    output:
        temp(directory("results/get-sra/{sample}-input"))
    shell:
        """
        wget -P {params.dir}/{wildcards.sample}-input {params.url}
        """

rule fastq_dump:
    input:
        rules.get_sra.output
    output:
        temp(directory("results/fastq-dump/{sample}"))
    params:
        run_id = lambda wc: basename(pep.get_sample(wc.sample).sra_accession[0])
    threads:
        16
    resources:
        time=60,
        mem=config.get("FASTERQDUMP_MEM", 8000),
        cpus=16
    log:
        "results/logs/fastq-dump/{sample}.log"
    conda:
        "../envs/sratools.yaml"
    shell:
       """
       fasterq-dump --mem {resources.mem}MB -s -S --include-technical -e {threads} -O {output} {input}/{params.run_id} 2> {log}
       """

rule fastq_dump_input:
    input:
        rules.get_sra_input.output
    output:
        temp(directory("results/fastq-dump/{sample}-input"))
    params:
        run_id = lambda wc: basename(pep.get_sample(wc.sample).sra_accession_input[0])
    threads:
        16
    resources:
        time=60,
        mem=config.get("FASTERQDUMP_MEM", 8000),
        cpus=16
    log:
        "results/logs/fastq-dump-input/{sample}.log"
    conda:
        "../envs/sratools.yaml"
    shell:
       """
       fasterq-dump --mem {resources.mem}MB -s -S --include-technical -e {threads} -O {output} {input}/{params.run_id} 2> {log}
       """


rule fastqc:
    input:
        rules.fastq_dump.output
    output:
        directory("results/fastqc/{sample}")
    conda:
        "../envs/fastqc.yaml"
    threads:
        16
    log:
        "results/logs/fastqc/{sample}.log"
    resources:
        time=20,
        cpus=16
    shell:
        """
        mkdir -p {output} &&
        fastqc -o {output} -t {threads} {input}/* 2> {log}
        """

rule bam_collate:
    input:
        "results/{file}.bam"
    output:
        "results/{file}.collate.bam"
    threads:
        16
    conda:
        "../envs/samtools.yaml"
    resources:
        time=20,
        cpus=16,
        mem=16000
    shell:
        """
        samtools collate -@ 16 -o {output} {input}
        """

rule bam_fixm:
    input:
        "results/{file}.collate.bam"
    output:
        "results/{file}.fixm.bam"
    conda:
        "../envs/samtools.yaml"
    threads:
        16
    resources:
        time=20,
        cpus=16,
        mem=16000
    shell:
        """
        samtools fixmate -@ {threads} -m {input} {output}
        """

rule bam_rmdup:
    input:
        "results/{file}.fixm.srt.bam"
    output:
        "results/{file}.rmdup.bam"
    conda:
        "../envs/samtools.yaml"
    threads:
        16
    resources:
        time=20,
        cpus=16,
        mem=16000
    shell:
        """
        samtools markdup -@ 16 -r {input} {output}
        """

rule bam_sort:
    input:
        "results/{file}.bam"
    output:
        "results/{file}.srt.bam"
    conda:
        "../envs/samtools.yaml"
    threads:
        16
    resources:
        time=20,
        cpus=16,
        mem=16000
    shell:
        """
        samtools sort -@ 16 {input} > {output}
        """

rule bam_index:
    input:
        "results/{file}.bam"
    output:
        "results/{file}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    threads:
        16
    resources:
        time=20,
        cpus=16,
        mem=16000
    shell:
        """
        samtools index -@ 16 {input}
        """
