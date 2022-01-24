#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2021, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import csv

SAMPLE_TABLE = "samples.tsv"
REF_GENOME_TABLE = "reference_genomes.tsv"
OUT_GENOME_TABLE = "outgroup_genomes.tsv"


def read_sample_list():
    """Read the user sample table"""

    samples = pd.read_csv(
        SAMPLE_TABLE,
        sep="\t",
        names=["Sample_Acc", "Species", "Cluster"],
    )

    return samples


def get_right_pathogen(wildcards, checkpoints):
    """Get the E. coli samples"""

    samples = read_sample_list()

    species = ""

    if wildcards.pathogen == "mdv":
        species = "MDV"

    patho_samples = samples[(samples["Species"] == species)]

    # check what cluster if necessary
    if hasattr(wildcards, "cluster") and (wildcards.cluster != "all"):
        patho_samples = patho_samples[patho_samples["Cluster"] == wildcards.cluster]

    inputs_all = []

    for key, sam in patho_samples.iterrows():
        inputs_all.append(sam["Sample_Acc"])

    return inputs_all


def get_ref_genome(wildcards):
    """Function to get the right ref genome for each fastbaps cluster"""

    ref_genomes = pd.read_csv(
        REF_GENOME_TABLE,
        sep="\t",
    )

    if not hasattr(wildcards, "cluster"):
        wildcards.cluster = "all"

    return ref_genomes.loc[
        ref_genomes["Cluster"] == wildcards.cluster, "Accession"
    ].to_list()[0]


def get_out_genome(wildcards):
    """Function to get the right ref genome for each fastbaps cluster"""

    ref_genomes = pd.read_csv(
        OUT_GENOME_TABLE,
        sep="\t",
    )

    if not hasattr(wildcards, "cluster"):
        wildcards.cluster = "all"

    return ref_genomes.loc[
        ref_genomes["Cluster"] == wildcards.cluster, "Accession"
    ].to_list()[0]


def get_genes(wildcards, checkpoints):
    """Get the paths to the correct individual gene alignments"""

    orfs = checkpoints.remove_orf_overlap.get(pathogen=wildcards.pathogen)
    bed = pd.read_csv(
        orfs.output[0],
        sep="\t",
        names=["chrom", "start", "end", "gene", "score", "strand"],
    )

    forward_genes = bed[bed["strand"] == "+"]
    reverse_genes = bed[bed["strand"] == "-"]

    gene_paths = []

    for key, gene in forward_genes.iterrows():
        gene_paths.append(
            f"genes_{wildcards.pathogen}/{gene['gene']}_{gene['start']}_{gene['end']}_aln.fasta"
        )

    for key, gene in reverse_genes.iterrows():
        gene_paths.append(
            f"genes_{wildcards.pathogen}/{gene['gene']}_{gene['start']}_{gene['end']}_aln_rev_comp.fasta"
        )

    return gene_paths


def get_intergenic(wildcards, checkpoints):
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

