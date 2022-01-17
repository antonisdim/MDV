#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"

MIN_FRAG_LEN = 0
MAX_FRAG_LEN = 1000

rule bowtie_index_accession:
    input:
        "refs/{accession}_mask.fasta",
    log:
        "refs/{accession}_mask_index.log",
    output:
        "refs/{accession}_mask.1.bt2l",
        "refs/{accession}_mask.2.bt2l",
        "refs/{accession}_mask.3.bt2l",
        "refs/{accession}_mask.4.bt2l",
        "refs/{accession}_mask.rev.1.bt2l",
        "refs/{accession}_mask.rev.2.bt2l",
    message:
        "Preparing the bowtie2 index for genome {wildcards.accession} of taxon MDV."
    threads: 4
    params:
        basename="refs/{accession}_mask",
    shell:
        "bowtie2-build --large-index --threads {threads} {input} {params.basename} &> {log}"


rule bowtie_aln_fasta:
    input:
        fasta="modern_genomes/{sample}_clean.split.500mer.5overlap.fasta",
        db_idx="refs/{accession}_mask.1.bt2l",
    log:
        "modern_genomes/{sample}_{accession}.log",
    output:
        bam_file="modern_genomes/{sample}_ref_{accession}.bam",
        bai_file="modern_genomes/{sample}_ref_{accession}.bam.bai",
    params:
        basename="refs/{accession}_mask",
    threads: 4
    message:
        "Aligning the reads from modern sample {wildcards.sample} against ref genome {wildcards.accession}."
    shell:
        "(bowtie2 -f --time --no-unal --no-mixed --threads {threads} -x {params.basename} -U {input.fasta} | "
        "samtools sort -O bam -o {output.bam_file} && samtools index {output.bam_file}) 2> {log}"


# todo add rule to align the raw modern data, maybe modify the samples.tsv file also with the help of two input functions, one for fasta one for fastq

rule bowtie_align_accession_paired_end:
    input:
        fastq_r1="modern_raw_reads/{sample}_R1_adRm.fastq.gz",
        fastq_r2="modern_raw_reads/{sample}_R2_adRm.fastq.gz",
        db_idx="refs/{accession}_mask.1.bt2l",
    log:
        "modern_raw_reads/{sample}_ref_{accession}.log",
    output:
        bam_file="modern_raw_reads/{sample}_ref_{accession}.bam",
        bai_file="modern_raw_reads/{sample}_ref_{accession}.bam.bai",
    params:
        basename="refs/{accession}_mask",
    threads: 4
    message:
        "Aligning the reads from sample {wildcards.sample} against ref genome {wildcards.accession}."
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "( bowtie2 --time --no-unal --no-mixed --threads {threads} "
        "-I {MIN_FRAG_LEN} -X {MAX_FRAG_LEN} -x {params.basename} "
        "-1 {input.fastq_r1} -2 {input.fastq_r2} "
        "| samtools sort -O bam -o {output.bam_file} && samtools index {output.bam_file} "
        ") 2> {log}"


