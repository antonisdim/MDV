#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"

import pandas as pd

from scripts.utilities import get_ref_genome, get_right_pathogen, get_out_genome


def get_modern_seqs_with_hvt(wildcards):
    """Get MLE sequences"""

    input_paths = []

    ref = get_ref_genome(wildcards)

    ncounts = checkpoints.count_fasta_n.get(pathogen=wildcards.pathogen)
    n_stats = pd.read_csv(ncounts.output[0], sep="\t", names=["sample", "percent"])
    n_stats_round = n_stats.round(0).sort_values(["percent"], ascending=False)
    selected_mdv = n_stats_round.loc[
        n_stats_round["percent"] <= 99.0,
    ]

    for key, sam in selected_mdv.iterrows():
        if not "OL" in sam["sample"]:
            input_paths.append(f"seqs_mdv/{sam['sample']}_ref_{ref}_CM.fasta")

    outgroup = get_out_genome(wildcards)
    input_paths.append(f"modern_genomes/{outgroup}.fasta")

    return input_paths


def get_modern_seqs_no_hvt(wildcards):
    """Get MLE sequences"""

    input_paths = []

    ref = get_ref_genome(wildcards)

    ncounts = checkpoints.count_fasta_n.get(pathogen=wildcards.pathogen)
    n_stats = pd.read_csv(ncounts.output[0], sep="\t", names=["sample", "percent"])
    n_stats_round = n_stats.round(0).sort_values(["percent"], ascending=False)
    selected_mdv = n_stats_round.loc[
        n_stats_round["percent"] <= 99.0,
    ]

    for key, sam in selected_mdv.iterrows():
        if not "OL" in sam["sample"] and sam["sample"] != "MG518371.1":
            input_paths.append(f"seqs_mdv/{sam['sample']}_ref_{ref}_CM.fasta")

    return input_paths


def get_ancient_seqs_beast(wildcards):
    """Get MLE sequences"""

    input_paths = []

    ref = get_ref_genome(wildcards)

    ncounts = checkpoints.count_fasta_n.get(pathogen=wildcards.pathogen)
    n_stats = pd.read_csv(ncounts.output[0], sep="\t", names=["sample", "percent"])
    n_stats_round = n_stats.round(0).sort_values(["percent"], ascending=False)
    selected_mdv = n_stats_round.loc[
        n_stats_round["percent"] <= 80.0,
    ]

    for key, sam in selected_mdv.iterrows():
        if "OL" in sam["sample"]:
            input_paths.append(f"seqs_mdv/{sam['sample']}_ref_{ref}_CM.fasta")

    return input_paths


def get_ancient_seqs_mle(wildcards):
    """Get MLE sequences"""

    input_paths = []
    ref = get_ref_genome(wildcards)

    ncounts = checkpoints.count_fasta_n.get(pathogen=wildcards.pathogen)
    n_stats = pd.read_csv(ncounts.output[0], sep="\t", names=["sample", "percent"])
    n_stats_round = n_stats.round(0).sort_values(["percent"], ascending=False)
    selected_mdv = n_stats_round.loc[
        (n_stats_round["percent"] <= 99.0) & (n_stats_round["percent"] > 80.0),
    ]

    for key, sam in selected_mdv.iterrows():
        input_paths.append(f"seqs_mdv/{sam['sample']}_ref_{ref}_CM.fasta")

    return input_paths


rule mafft_modern:
    input:
        get_modern_seqs_with_hvt,
    log:
        "aln_{pathogen}/mdv_modern_plus_HVT_aln.log",
    output:
        non_aln=temp("aln_{pathogen}/mdv_modern_plus_HVT.fasta"),
        aln="aln_{pathogen}/mdv_modern_plus_HVT_aln.fasta",
    message:
        "Aligning the modern MDV and HVT seqs with mafft."
    threads: workflow.cores
    shell:
        "cat {input} > {output.non_aln} &&"
        "mafft --maxiterate 1000 --thread {threads} --nwildcard {output.non_aln} 1> {output.aln} 2> {log} && "
        "sed -i 's/ Meleagrid herpesvirus 1, complete genome//g' {output.aln}"


rule mafft_beast:
    input:
        mod_aln="aln_{pathogen}/mdv_modern_plus_HVT_aln.fasta",
        beast_anc=get_ancient_seqs_beast,
    log:
        "aln_{pathogen}/mdv_mod_anc_plus_HVT_aln_BEAST.log",
    output:
        "aln_{pathogen}/mdv_mod_anc_plus_HVT_aln_BEAST.fasta",
    message:
        "Adding the ancient genomes for the BEAST MLE tree alignment."
    threads: workflow.cores
    shell:
        "(cp {input.mod_aln} aln_{wildcards.pathogen}/mdv_mod_anc_plus_HVT_aln_tmp_0.fasta && "
        "count=0 && "
        "for seq in {input.beast_anc}; "
        "do res=$((count+1)); "
        "mafft --maxiterate 1000 --thread {threads} --nwildcard "
        "--add $seq aln_{wildcards.pathogen}/mdv_mod_anc_plus_HVT_aln_tmp_${{count}}.fasta  1> "
        "aln_{wildcards.pathogen}/mdv_mod_anc_plus_HVT_aln_tmp_${{res}}.fasta; "
        "count=$((count+1)); done && "
        "mv aln_{wildcards.pathogen}/mdv_mod_anc_plus_HVT_aln_tmp_${{count}}.fasta {output} && "
        "for tmp in aln_{wildcards.pathogen}/mdv_mod_anc_plus_HVT_aln_tmp_*.fasta; do unlink $tmp; done) 2> {log}"


rule mafft_mle:
    input:
        beast_aln="aln_{pathogen}/mdv_mod_anc_plus_HVT_aln_BEAST.fasta",
        mle_anc=get_ancient_seqs_mle,
    log:
        "aln_{pathogen}/mdv_mod_anc_plus_HVT_aln_MLE.log",
    output:
        "aln_{pathogen}/mdv_mod_anc_plus_HVT_aln_MLE.fasta",
    message:
        "Adding the ancient genomes for the RAxML MLE tree alignment."
    threads: workflow.cores
    shell:
        "(cp {input.beast_aln} aln_{wildcards.pathogen}/mdv_mod_anc_plus_HVT_aln_tmp_0.fasta && "
        "count=0 && "
        "for seq in {input.mle_anc}; "
        "do res=$((count+1)); "
        "mafft --maxiterate 1000 --thread {threads} --nwildcard "
        "--add $seq aln_{wildcards.pathogen}/mdv_mod_anc_plus_HVT_aln_tmp_${{count}}.fasta  1> "
        "aln_{wildcards.pathogen}/mdv_mod_anc_plus_HVT_aln_tmp_${{res}}.fasta; "
        "count=$((count+1)); done && "
        "mv aln_{wildcards.pathogen}/mdv_mod_anc_plus_HVT_aln_tmp_${{count}}.fasta {output} && "
        "for tmp in aln_{wildcards.pathogen}/mdv_mod_anc_plus_HVT_aln_tmp_*.fasta; do unlink $tmp; done) 2> {log}"


rule beast_aln:
    input:
        modern=get_modern_seqs_no_hvt,
        ancient=get_ancient_seqs_beast,
    output:
        "aln_{pathogen}/mdv_mod_anc_no_HVT_aln_BEAST.fasta",
    message:
        "Creating multifasta alignment file for BEAST."
    shell:
        "cat {input.modern} {input.ancient} > {output}"
