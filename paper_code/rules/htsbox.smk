#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"


rule consensus_fasta_ancient:
    input:
        bam_file="merged_bams/{sample}_adRm_sorted_rmdup_m.bam",
        ref="refs/{accession}_mask.fasta",
        faidx="refs/{accession}_mask.fasta.fai",
    log:
        "seqs_{pathogen}/{sample}_ref_{accession}.log",
    output:
        "seqs_{pathogen}/{sample}_ref_{accession}_CM.fasta",
    message:
        "Creating the consensus fasta sequence for sample {wildcards.sample} for "
        "{wildcards.pathogen}."
    shell:
        "(htsbox pileup -f {input.ref} -l 25 -T 5 -q 30 -Q 30 -M -s 5 {input.bam_file} 1> {output} && "
        "sed -i 's/>EF523390.1/>{wildcards.sample}/g' {output}) 2> {log}"


rule consensus_fasta_modern:
    input:
        bam_file="modern_genomes/{sample}_ref_{accession}.bam",
        ref="refs/{accession}_mask.fasta",
        faidx="refs/{accession}_mask.fasta.fai",
    log:
        "seqs_{pathogen}/{sample}_ref_{accession}.log",
    output:
        "seqs_{pathogen}/{sample}_ref_{accession}_CM.fasta",
    message:
        "Creating the consensus fasta sequence for sample {wildcards.sample} for "
        "{wildcards.pathogen}."
    shell:
        "(htsbox pileup -f {input.ref} -l 25 -M {input.bam_file} 1> {output} && "
        "sed -i 's/>EF523390.1/>{wildcards.sample}/g' {output}) 2> {log}"


rule consensus_fasta_modern_raw:
    input:
        bam_file="modern_raw_reads/{sample}_ref_{accession}.bam",
        ref="refs/{accession}_mask.fasta",
        faidx="refs/{accession}_mask.fasta.fai",
    log:
        "seqs_{pathogen}/{sample}_ref_{accession}.log",
    output:
        "seqs_{pathogen}/{sample}_ref_{accession}_CM.fasta",
    message:
        "Creating the consensus fasta sequence for sample {wildcards.sample} for "
        "{wildcards.pathogen}."
    shell:
        "(htsbox pileup -f {input.ref} -l 25 -T 3 -q 30 -Q 30 -M -s 5 {input.bam_file} 1> {output} && "
        "sed -i 's/>EF523390.1/>{wildcards.sample}/g' {output}) 2> {log}"
