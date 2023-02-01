#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_ref_genome, get_right_pathogen

import pandas as pd

checkpoint remove_orf_overlap:
    input:
        "aux_files/{pathogen}_gene_loci.bed",
    output:
        "aux_files/{pathogen}_gene_loci_no_overlap.bed",
    message:
        "Removing overlapping ORF regions from the MDV genome."
    script:
        "../scripts/remove_orf_overlap.py"


def get_ref_fai(wildcards):
    """Get the correct fasta index"""

    ref = get_ref_genome(wildcards)

    return f"refs/{ref}_mask.fasta.fai"


rule intergenic_bed:
    input:
        bed="aux_files/{pathogen}_gene_loci.bed",
        fai=get_ref_fai,
    output:
        temp_bed=temp("aux_files/{pathogen}_gene_merged_loci_tmp.bed"),
        gen_file=temp("aux_files/{pathogen}_mask_gen_file.txt"),
        bed=temp("aux_files/{pathogen}_gene_non_coding.bed"),
    message:
        "Finding the coordinates of the intergenic regions of the {wildcards.pathogen} genome."
    shell:
        "bedtools merge -i {input.bed} > {output.temp_bed} && "
        "cat {input.fai} | cut -f 1-2 > {output.gen_file} && "
        "bedtools complement -i {output.temp_bed} -g {output.gen_file} > {output.bed}"


checkpoint fix_bed_coord:
    input:
        "aux_files/{pathogen}_gene_non_coding.bed",
    output:
        "aux_files/{pathogen}_gene_non_coding_1idx.bed",
    message:
        "Getting the intergenic regions bed file with 1-indexed coordinates for {wildcards.pathogen}."
    script:
        "../scripts/fix_bed_coord.py"


rule indiv_genes:
    input:
        "aln_{pathogen}/mdv_mod_anc_no_HVT_aln_BEAST.fasta",
    output:
        "genes_{pathogen}/{gene}_{start}_{end}_aln.fasta",
    message:
        "Exatracting the alignments for genomic region {wildcards.gene} "
        "from {wildcards.start} to {wildcards.end}."
    shell:
        "seqkit subseq -r {wildcards.start}:{wildcards.end} {input} > {output}"


rule reverse_comp:
    input:
        "genes_{pathogen}/{gene}_{start}_{end}_aln.fasta",
    output:
        "genes_{pathogen}/{gene}_{start}_{end}_aln_rev_comp.fasta",
    message:
        "Fixing the directionality of ORF {wildcards.gene}."
    shell:
        "seqtk seq -r -l 70 {input} > {output}"


def get_beast_regions(wildcards):
    """Get the paths to the correct genomic regions alignments for BEAST"""

    alns = checkpoints.percent_overlap_orfs.get(
        pathogen=wildcards.pathogen, region=wildcards.region
    )
    paths = pd.read_csv(alns.output[0], sep="\t", names=["path"])

    region_paths = []

    for key, gene in paths.iterrows():
        region_paths.append(gene["path"])

    return region_paths


rule beast_regions:
    input:
        get_beast_regions,
    output:
        "aln_{pathogen}/{pathogen}_{region}_region_aln.fasta",
    message:
        "Creating the alignment for the {wildcards.region} regions of {wildcards.pathogen}."
    shell:
        "seqkit concat {input} > {output}"


rule extract_codon:
    input:
        "aln_{pathogen}/{pathogen}_{region}_region_aln.fasta",
    output:
        "aln_{pathogen}/{pathogen}_{region}_region_codon{codon}_aln.fasta",
    message:
        "Extracting codon position {wildcards.codon} for the {wildcards.pathogen} BEAST "
        "{wildcards.region} aln."
    script:
        "../scripts/extract_codon.py"
