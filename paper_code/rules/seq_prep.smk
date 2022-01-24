#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_ref_genome, get_right_pathogen

import pandas as pd


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


def get_genes(wildcards):
    """Get the paths to the correct individual gene alignments"""

    gene_paths = []

    if wildcards.region == "coding":

        orfs = checkpoints.remove_orf_overlap.get(pathogen=wildcards.pathogen)
        bed = pd.read_csv(
            orfs.output[0],
            sep="\t",
            names=["chrom", "start", "end", "gene", "score", "strand"],
        )

        forward_genes = bed[bed["strand"] == "+"]
        reverse_genes = bed[bed["strand"] == "-"]

        for key, gene in forward_genes.iterrows():
            gene_paths.append(
                f"genes_{wildcards.pathogen}/{gene['gene']}_{gene['start']}_{gene['end']}_aln.fasta"
            )

        for key, gene in reverse_genes.iterrows():
            gene_paths.append(
                f"genes_{wildcards.pathogen}/{gene['gene']}_{gene['start']}_{gene['end']}_aln_rev_comp.fasta"
            )

    elif wildcards.region == "noncoding":

        orfs = checkpoints.fix_bed_coord.get(pathogen=wildcards.pathogen)
        bed = pd.read_csv(
            orfs.output[0],
            sep="\t",
            names=["chrom", "start", "end"],
        )

        for key, gene in bed.iterrows():
            gene_paths.append(
                f"genes_{wildcards.pathogen}/{gene['chrom']}_{gene['start']}_{gene['end']}_aln.fasta"
            )

    return gene_paths


checkpoint percent_overlap_orfs:
    input:
        get_genes,
    output:
        "aux_files/{pathogen}_{region}_region.txt",
    message:
        "Selecting {wildcards.region} regions with enough sequence overlap in their alignments."
    script:
        "../scripts/percent_overlap_orfs.py"


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


def get_intergenic(wildcards):
    """Get the paths to the correct individual gene alignments"""

    orfs = checkpoints.fix_bed_coord.get(pathogen=wildcards.pathogen)
    bed = pd.read_csv(
        orfs.output[0],
        sep="\t",
        names=["chrom", "start", "end"],
    )

    gene_paths = []

    for key, gene in bed.iterrows():
        gene_paths.append(
            f"genes_{wildcards.pathogen}/{gene['gene']}_{gene['start']}_{gene['end']}_aln.fasta"
        )

    return gene_paths
