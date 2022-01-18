#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


rule adapterremoval_paired_end:
    input:
        fastq_r1="modern_raw_reads/{accession}_R1.fastq.gz",
        fastq_r2="modern_raw_reads/{accession}_R2.fastq.gz",
    log:
        "modern_raw_reads/{accession}_adRm.log",
    output:
        fastq_r1="modern_raw_reads/{accession}_R1_adRm.fastq.gz",
        fastq_r2="modern_raw_reads/{accession}_R2_adRm.fastq.gz",
    message:
        "Trimming sequencing adapters from files {input.fastq_r1} and {input.fastq_r2}."
    threads: 2
    params:
        basename="modern_raw_reads/{accession}",
    shell:
        "(AdapterRemoval"
        "   --file1 {input.fastq_r1}"
        "   --file2 {input.fastq_r2} "
        "   --basename {params.basename}"
        "   --gzip "
        "   --minlength 15 "
        "   --threads {threads}"
        "   --trimns &&"
        " mv {params.basename}.pair1.truncated.gz {output.fastq_r1} && "
        " mv {params.basename}.pair2.truncated.gz {output.fastq_r2} "
        ") 2> {log}"
