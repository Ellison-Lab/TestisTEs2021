rule get_sra_10x:
    params:
        url = lambda wc: pep.get_sample(wc.sample).sra_accession,
        dir = "results/get-sra-10x/"
    output:
        temp(directory("results/get-sra-10x/{sample}"))
    shell:
        """
        wget -P {params.dir} {params.url}
        """

checkpoint fastq_dump_10x:
    input:
        rules.get_sra_10x.output
    output:
        temp(directory("results/fastq-dump-10x/{sample}/"))
    params:
        run_id = lambda wc: basename(pep.get_sample(wc.sample).sra_accession[0])
    threads:
        4
    resources:
        time=60,
        mem=8000,
        cpus=4
    log:
        "results/logs/fastq-dump-10x/{sample}.log"
    conda:
        "../envs/sratools.yaml"
    shell:
       """
       fasterq-dump -s -S --include-technical -e {threads} -O {output} {input}/{params.run_id} 2> {log}
       """

def determine_layout(wildcards, dir = "results/fastq-dump-10x/"):
    checkpoint_output = checkpoints.fastq_dump_10x.get(**wildcards).output[0]
    run_id = basename(pep.get_sample(wildcards.sample).sra_accession)
    base = "{d}/{s}/{r}".format(d=dir,s=wildcards.sample, r=run_id)
    glc = glob_wildcards(base + "_{i}.fastq")
    op = expand("{b}_{i}.fastq", b=base,i=glc.i)
    return(op)

rule rename_fqs_10x:
    input:
        determine_layout
    output:
        "results/fastq-rename-10x/{sample}_S1_R1_001.fastq",
        "results/fastq-rename-10x/{sample}_S1_R2_001.fastq",
    log:
        "logs/fastq-rename-10x/{sample}.log"
    run:
        s = wildcards.sample
        run_id = basename(pep.get_sample(wc.sample).sra_accession)
        if (len(input) == 2):
            shell("cp {s} {o}".format(s=input[0], o=output[0])),
            shell("cp {s} {o}".format(s=input[1], o=output[1]))
        elif (len(input) == 3):
            shell("cp {s} {o}".format(s=input[1], o=output[0])),
            shell("cp {s} {o}".format(s=input[2], o=output[1]))

rule cellranger_mkref:
    input:
        fa = rules.combine_fas.output,
        gtf = rules.fix_strand_gtf.output
    output:
        directory("results/cellranger/idx")
    resources:
        time="18:00:00",
        mem=int(config.get("CELLRANGER_MEM",16000)/1000.0),
        cpus=1
    shell:
        """
        ~/cellranger-3.1.0/cellranger mkref --genome={output} \
            --fasta={input.fa} \
	    --genes={input.gtf} \
	    --memgb=16
        """

rule cellranger_count:
    input:
        idx = rules.cellranger_mkref.output,
        fqs = rules.rename_fqs_10x.output,
    output:
        directory("results/cellranger/counts/{sample}")
    threads:
        config.get("CELLRANGER_COUNT_CPUS",32)
    params:
        ver = config.get("CELLRANGER_VERSION",'3.1.0')
    resources:
        time="18:00:00",
        mem=int(config.get("CELLRANGER_MEM",128000)/1000.0),
        cpus=config.get("CELLRANGER_CPUS",32),
    shell:
        """
        ~/cellranger-{params.ver}/cellranger count --id={output} \
            --sample={wildcards.sample} \
            --fastqs=results/fastq-rename-10x/{wildcards.sample} \
	    --transcriptome={input.idx} \
	    --localcores={threads} \
	    --localmem={resources.mem}
        """
