#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"

import os


rule run_raxml_gtr_gamma:
    input:
        "aln_{pathogen}/mdv_mod_anc_{outgroup}_HVT_aln_{cluster}.fasta",
    log:
        "trees_{pathogen}/RAxML_{pathogen}_{outgroup}_HVT_{cluster}.log",
    output:
        best_tree="trees_{pathogen}/RAxML_bestTree.{pathogen}_{outgroup}_HVT_{cluster}",
        bipartition_labels="trees_{pathogen}/RAxML_bipartitionsBranchLabels.{pathogen}_{outgroup}_HVT_{cluster}",
        bipartition=(
            "trees_{pathogen}/RAxML_bipartitions.{pathogen}_{outgroup}_HVT_{cluster}"
        ),
        bootstrap="trees_{pathogen}/RAxML_bootstrap.{pathogen}_{outgroup}_HVT_{cluster}",
        info="trees_{pathogen}/RAxML_info.{pathogen}_{outgroup}_HVT_{cluster}",
    message:
        "Running raxml for {wildcards.pathogen} ({wildcards.outgroup} HVT) for {wildcards.cluster}."
    threads: workflow.cores
    params:
        basename="{pathogen}_{outgroup}_HVT_{cluster}",
        workdir=lambda wildcards: os.path.abspath(f"trees_{wildcards.pathogen}"),
    shell:
        "(raxmlHPC-PTHREADS -f a -x 12345 -p 12345 -T {threads} -m GTRGAMMA -k -# 100 "
        "-s {input} -n {params.basename} -w {params.workdir}) 2> {log}"


rule prune_tree:
    input:
        "trees_{pathogen}/RAxML_bipartitions.{pathogen}_plus_HVT_{cluster}",
    output:
        "trees_{pathogen}/RAxML_bipartitions.{pathogen}_pruned_HVT_{cluster}",
    message:
        "Pruning the HVT root from the {wildcards.pathogen} tree for {wildcards.cluster}."
    shell:
        "Rscript scripts/prune_tree.R {input} {output}"
