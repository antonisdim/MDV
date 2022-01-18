#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_ref_genome, get_right_pathogen


rule clean_from_BACs:
    input:
        "modern_genomes/{sample}.fasta",
    output:
        temp("modern_genomes/{sample}_clean.fasta"),
    message:
        "Cleaning up any residual BAC sequences in sample {wildcards.sample}"
    script:
        "../scripts/clean_from_BACs.py"


rule pyfasta_split:
    input:
        "modern_genomes/{sample}_clean.fasta",
    output:
        fasta=temp("modern_genomes/{sample}_clean.split.500mer.5overlap.fasta"),
        flat=temp("modern_genomes/{sample}_clean.fasta.flat"),
        gdx=temp("modern_genomes/{sample}_clean.fasta.gdx"),
    message:
        "Fragmenting sample {wildcards.sample} into kmers of 500 bp with a 5 bp overalp."
    shell:
        "pyfasta split -n 1 -k 500 -o 5 {input}"


def get_clean_masked_alns(wildcards):
    """Get the samples belonging to a cluster for alignments"""

    samples_list = get_right_pathogen(wildcards, checkpoints)
    input_paths = []

    ref = get_ref_genome(wildcards)

    for sample in samples_list:
        input_paths.append(f"seqs_mdv/{sample}_ref_{ref}_CM.fasta")

    return input_paths


rule concat_alns:
    input:
        get_clean_masked_alns,
    output:
        temp("aln_{pathogen}/all_seqs_plus_HVT_CM.fasta"),
    message:
        "Concatenating all the clean masked aligned samples into a mutlifasta files."
    shell:
        "cat {input} > {output}"


checkpoint count_fasta_n:
    input:
        "aln_{pathogen}/all_seqs_plus_HVT_CM.fasta",
    output:
        "aln_{pathogen}/all_seqs_plus_HVT_CM_nstats.tsv",
    message:
        "Calculating the percentage of Ns in the consensus clean and masked fasta sequences."
    shell:
        "seqtk comp {input} | column -t | awk '{{print $1\"\t\"$9/$2*100}}' > {output}"
