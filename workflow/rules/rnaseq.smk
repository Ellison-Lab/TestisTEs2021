rule cat_fqs:
    input:
        rules.fastq_dump.output
    output:
        r1=temp("results/rnaseq/fastq-combine/{sample}/{sample}_r1.concat.fastq.gz"),
        r2=temp("results/rnaseq/fastq-combine/{sample}/{sample}_r2.concat.fastq.gz")
    params:
        r1 = lambda wc: ["results/fastq-dump/{s}/".format(s=wc.sample) + basename(x) + "_1.fastq" for x in pep.get_sample(wc.sample).sra_accession],
        r2 = lambda wc: ["results/fastq-dump/{s}/".format(s=wc.sample) + basename(x) + "_2.fastq" for x in pep.get_sample(wc.sample).sra_accession],
    shell:
        """
        cat {params.r1} > {output.r1} && \
        cat {params.r2} > {output.r2}
        """

rule trim_pe:
    input:
        fq1 = rules.cat_fqs.output.r1,
        fq2 = rules.cat_fqs.output.r2
    output:
        r1 = "results/trimmed/{sample}/{sample}_r1.trimmed.fq.gz",
        r2 = "results/trimmed/{sample}/{sample}_r2.trimmed.fq.gz",
        html = "results/trimmed/{sample}/{sample}_fastp.html",
        json = "results/trimmed/{sample}/{sample}_fastp.json"
    threads:
        2
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp --in1 {input.fq1} --in2 {input.fq2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.sample}_fastp"

rule generate_star_idx:
    """
    no dedup here: https://www.nature.com/articles/srep25533
    https://kb.10xgenomics.com/hc/en-us/articles/115004415786-What-parameters-are-used-for-STAR-alignment-
    https://github.com/10XGenomics/cellranger/blob/master/mro/stages/counter/align_reads/__init__.py
    https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/reference.py
    """
    input:
        genome = "results/custom-genome/combined.fasta",
        gtf = "results/custom-genome/combined.fixed.gtf"
    output:
        "results/rnaseq/index/chrLength.txt",
        "results/rnaseq/index/chrNameLength.txt",
        "results/rnaseq/index/chrName.txt",
        "results/rnaseq/index/chrStart.txt",
        "results/rnaseq/index/exonGeTrInfo.tab",
        "results/rnaseq/index/exonInfo.tab",
        "results/rnaseq/index/Genome",
        "results/rnaseq/index/geneInfo.tab",
        "results/rnaseq/index/genomeParameters.txt",
        "results/rnaseq/index/SA",
        "results/rnaseq/index/SAindex",
        "results/rnaseq/index/sjdbList.fromGTF.out.tab",
        "results/rnaseq/index/sjdbInfo.txt",
        "results/rnaseq/index/sjdbList.out.tab",
        "results/rnaseq/index/transcriptInfo.tab",
    resources:
        time=40,
        mem=24000,
        cpus=4
    threads:
        4
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode genomeGenerate "
        "--genomeDir results/rnaseq/index/ --runThreadN {threads} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbGTFfile {input.gtf} "
        "--genomeSAindexNbases 11"

rule rnaseq_star_align:
    input:
        rules.generate_star_idx.output,
        reads1 = rules.trim_pe.output.r1,
        reads2 = rules.trim_pe.output.r2,
        gtf = "results/custom-genome/combined.fixed.gtf"
    output:
        directory("results/rnaseq/star/{sample}/")
    threads:
        12
    resources:
        time=120,
        mem=24000,
        cpus=12
    conda:
        "../envs/star.yaml"
    shell:
        """
        mkdir -p {output}
        STAR \
        --genomeDir results/rnaseq/index/ --runThreadN {threads} \
        --readFilesIn {input.reads1} {input.reads2} \
        --outSAMmultNmax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --sjdbGTFfile {input.gtf} \
        --readFilesCommand zcat \
        --outFileNamePrefix {output}/ \
        --chimSegmentMin 12 \
        --twopassMode Basic \
          --chimJunctionOverhangMin 8 \
          --chimOutJunctionFormat 1 \
          --alignSJDBoverhangMin 10 \
          --alignMatesGapMax 100000 \
          --alignIntronMax 100000 \
          --alignSJstitchMismatchNmax 5 -1 5 5 \
          --outSAMattrRGline ID:GRPundef \
          --chimMultimapScoreRange 3 \
          --chimScoreJunctionNonGTAG 0 \
          --chimMultimapNmax 20 \
          --chimNonchimScoreDropMin 10 \
          --peOverlapNbasesMin 12 \
          --peOverlapMMp 0.1 \
          --alignInsertionFlush Right \
          --alignSplicedMateMapLminOverLmate 0 \
          --alignSplicedMateMapLmin 30 \
        --chimOutType Junctions
        """
#
#
# rule index_bam:
#     input:
#         "results/star/{sample}/"
#     output:
#         touch('results/bam.index.done')
#     threads:
#         1
#     conda:
#         "envs/bowtie2.yaml"
#     shell:
#         "samtools index {input}/Aligned.sortedByCoord.out.bam"
