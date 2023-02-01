#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"


from scripts.utilities import get_ref_genome, read_mdv_seq_stats

EXCLUSION_LIST = [
    "JQ806362.1",
    "JQ809691.1",
    "JQ809692.1",
    "JQ820250.1",
    "JQ836662.1",
    "KT833851.1",
    "KT833852.1",
    "KX290015.1",
    "KX290016.1",
    "NC_002229.3",
    "AF147806.2",
    "MG518371.1",
]


def get_modern_seqs(wildcards):
    """Get MLE sequences"""

    input_paths = []
    ref = get_ref_genome(wildcards)
    n_stats_round = read_mdv_seq_stats(wildcards, checkpoints)

    selected_mdv = n_stats_round.loc[
        n_stats_round["percent"] <= 99.0,
    ]

    if wildcards.tree == "BEAST":
        selected_mdv = selected_mdv[~selected_mdv["sample"].isin(EXCLUSION_LIST)]

    for key, sam in selected_mdv.iterrows():
        if not "OL" in sam["sample"]:
            input_paths.append(f"seqs_mdv/{sam['sample']}_ref_{ref}_CM.fasta")

    return input_paths


def get_ancient_seqs(wildcards):
    """Get MLE sequences"""

    input_paths = []
    ref = get_ref_genome(wildcards)
    n_stats_round = read_mdv_seq_stats(wildcards, checkpoints)

    if wildcards.tree == "BEAST":
        selected_mdv = n_stats_round.loc[
            n_stats_round["percent"] <= 80.0,
        ]
    else:
        selected_mdv = n_stats_round.loc[
            n_stats_round["percent"] <= 99.0,
        ]

    for key, sam in selected_mdv.iterrows():
        if "OL" in sam["sample"]:
            input_paths.append(f"seqs_mdv/{sam['sample']}_ref_{ref}_CM.fasta")

    return input_paths


rule tree_alignments:
    input:
        modern=get_modern_seqs,
        ancient=get_ancient_seqs,
    output:
        "aln_{pathogen}/mdv_mod_anc_{outgroup}_HVT_aln_{tree}.fasta",
    message:
        "Creating multifasta alignment file for {wildcards.tree}."
    shell:
        "cat {input.modern} {input.ancient} > {output}"
