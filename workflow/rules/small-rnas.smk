
"""

https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
"""

rule shortstack_ebwt:
    input:
        ref =rules.combine_fas.output
    output:
        directory('results/small-rna/ebwt/')
    threads:
        12
    conda:
        '../envs/shortstack.yaml'
    resources:
        mem=config.get('SHORTSTACK_MEM'),
        cpus=12,
        time=20,
    shell:
        """
        mkdir -p {output} &&
        cp {input.ref} {output} &&
        bowtie-build --threads {threads} {input.ref} {output}/combined
        """

rule shortstack:
    input:
        fq = rules.fastq_dump.output,
        idx = rules.shortstack_ebwt.output
    output:
        directory('results/small-rna/shortstack/{sample}')
    params:
        dir = lambda wc: "results/small-rna/shortstack/{s}/".format(s=wc.sample),
        reads = lambda wc: "results/fastq-dump/{s}/{r}.fastq".format(s=wc.sample, r= basename(pep.get_sample(wc.sample).sra_accession[0])),
        mem = int(config.get('SHORTSTACK_MEM')/1000),
        adap = lambda wc: "--adapter {a}".format(a=pep.get_sample(wc.sample).adapter[0]) if (pep.get_sample(wc.sample).adapter[0] != "NA") else ""
    threads:
        20
    conda:
        '../envs/shortstack.yaml'
    resources:
        mem=config.get('SHORTSTACK_MEM'),
        cpus=22,
        time=60,
    shell:
        """
        ShortStack --readfile {params.reads} --genomefile \
            {input.idx}/combined.fasta --outdir {output} --sort_mem {params.mem}G \
            --bowtie_cores {threads} \
            {params.adap}
        """

rule cleanup_shortstack:
    input:
        rules.shortstack.output
    output:
        "results/small-rna/shortstack-clean/{sample}/{sample}.bam"
    threads:
        22
    conda:
        '../envs/shortstack.yaml'
    resources:
        mem=config.get('SHORTSTACK_MEM'),
        cpus=22,
        time=20,
    params:
        alns = lambda wc: "results/small-rna/shortstack/{s}/{r}.bam".format(s=wc.sample,r= basename(pep.get_sample(wc.sample).sra_accession[0])),
    shell:
        """
        samtools view -O BAM {params.alns} > {output}
        """
