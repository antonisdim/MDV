#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"


import pandas as pd


checkpoint count_fasta_n:
    input:
        "aln_{pathogen}/all_seqs_plus_HVT_CM.fasta",
    output:
        "aln_{pathogen}/all_seqs_plus_HVT_CM_nstats.tsv",
    message:
        "Calculating the percentage of Ns in the consensus clean and masked fasta sequences."
    shell:
        "seqtk comp {input} | column -t | awk '{{print $1\"\t\"$9/$2*100}}' > {output}"


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
