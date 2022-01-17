#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_ref_genome, get_right_pathogen


rule clean_from_BACs:
    input:
        "modern_genomes/{sample}.fasta"
    output:
        temp("modern_genomes/{sample}_clean.fasta")
    message:
        "Cleaning up any residual BAC sequences in sample {wildcards.sample}"
    script:
        "../scripts/clean_from_BACs.py"


rule pyfasta_split:
    input:
        "modern_genomes/{sample}_clean.fasta"
    output:
        fasta=temp("modern_genomes/{sample}_clean.split.500mer.5overlap.fasta"),
        flat=temp("modern_genomes/{sample}_clean.fasta.flat"),
        gdx=temp("modern_genomes/{sample}_clean.fasta.gdx"),
    message:
        "Fragmenting sample {wildcards.sample} into kmers of 500 bp with a 5 bp overalp."
    shell:
        "pyfasta split -n 1 -k 500 -o 5 {input}"


# todo do the htsbox for OL samples too
