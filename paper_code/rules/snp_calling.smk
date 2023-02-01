#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos"
__copyright__ = "Copyright 2022, University of Oxford"
__email__ = "ea.dimopoulos@gmail.com"
__license__ = "MIT"

from scripts.utilities import get_ref_genome, get_right_pathogen


rule ref_dict:
    input:
        "refs/{accession}_mask.fasta",
    output:
        dict="refs/{accession}_mask.dict",
        fai="refs/{accession}_mask.fasta.fai",
    message:
        "Preparing the sequene dict for {wildcards.accession}."
    shell:
        "picard CreateSequenceDictionary --REFERENCE {input} --OUTPUT {output.dict} && "
        "samtools faidx {input}"


rule trim_bams:
    input:
        "merged_bams/{sample}_adRm_sorted_rmdup_m.bam",
    output:
        temp("merged_bams/{sample}_adRm_sorted_rmdup_m_trim.bam"),
    message:
        "Trimming 5 bp from each of the reads from sample {wildcards.sample} "
        "to prepare them for SNP calling."
    shell:
        "bam trimBam {input} {output} -L 5 -R 5"


rule fix_read_groups:
    input:
        "merged_bams/{sample}_adRm_sorted_rmdup_m_trim.bam",
    log:
        "trimmed_bams_snp/{sample}_fix_read_groups.log",
    output:
        bam=temp("trimmed_bams_snp/{sample}_adRm_sorted_rmdup_m_trim_RG.bam"),
        bai=temp("trimmed_bams_snp/{sample}_adRm_sorted_rmdup_m_trim_RG.bam.bai"),
    message:
        "Fixing the read group headers in sample {wildcards.sample}."
    shell:
        "(picard AddOrReplaceReadGroups --INPUT {input} --OUTPUT {output.bam} "
        "--RGID M_A00363:131 --RGLB lib1 --RGPL Illumina --RGPU HHY2CDSXX:2 "
        "--RGSM {wildcards.sample} --VALIDATION_STRINGENCY LENIENT && "
        "samtools index {output.bam}) 2> {log}"


rule haplotype_caller:
    input:
        bam_file="trimmed_bams_snp/{sample}_adRm_sorted_rmdup_m_trim_RG.bam",
        bai_idx="trimmed_bams_snp/{sample}_adRm_sorted_rmdup_m_trim_RG.bam.bai",
        ref="refs/{accession}_mask.fasta",
        ref_dict="refs/{accession}_mask.dict",
    log:
        "gatk_mdv/{sample}_ref_{accession}.log",
    output:
        gvcf="gatk_mdv/{sample}_ref_{accession}.g.vcf.gz",
        tbi="gatk_mdv/{sample}_ref_{accession}.g.vcf.gz.tbi",
    message:
        "Preparing the GVCF file for MDV sample {wildcards.sample}."
    threads: 1
    shell:
        "(gatk HaplotypeCaller --input {input.bam_file} --output {output.gvcf} "
        "--reference {input.ref} --min-base-quality-score 30 "
        "--output-mode EMIT_ALL_ACTIVE_SITES --sample-ploidy 1 -ERC GVCF "
        "--do-not-run-physical-phasing true) 2> {log}"


def get_cluster_gvcfs(wildcards):
    """Get the samples belonging to a cluster for alignments"""

    samples_list = get_right_pathogen(wildcards)
    input_paths = []

    ref = get_ref_genome(wildcards)

    for sample in samples_list:
        if "OL" in sample:
            input_paths.append(f"gatk_{wildcards.pathogen}/{sample}_ref_{ref}.g.vcf.gz")

    return input_paths


def get_ref_fasta(wildcards):
    """Get the correct fasta index"""

    ref = get_ref_genome(wildcards)

    return f"refs/{ref}_mask.fasta"


rule combine_gvcfs:
    input:
        gvcfs=get_cluster_gvcfs,
        ref=get_ref_fasta,
    log:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_var_cohort.log",
    output:
        gvcf="gatk_{pathogen}/{pathogen}_cluster_{cluster}_var_cohort.g.vcf.gz",
    message:
        "Combining GVCF files for {wildcards.pathogen} cluster {wildcards.cluster}."
    wrapper:
        "0.80.2/bio/gatk/combinegvcfs"


rule genotype_gvcfs:
    input:
        gvcf="gatk_{pathogen}/{pathogen}_cluster_{cluster}_var_cohort.g.vcf.gz",
        ref=get_ref_fasta,
    log:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_var_joint_raw.log",
    output:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_var_joint_raw.vcf.gz",
    message:
        "Genotyping the cohort of GVCF files for {wildcards.pathogen} cluster {wildcards.cluster}."
    shell:
        "(gatk GenotypeGVCFs --reference {input.ref} --variant {input.gvcf} --output {output}) 2> {log}"


rule select_snps:
    input:
        vcf="gatk_{pathogen}/{pathogen}_cluster_{cluster}_var_joint_raw.vcf.gz",
        ref=get_ref_fasta,
    log:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_joint_raw.log",
    output:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_joint_raw.vcf.gz",
    message:
        "Slecting only the SNP variants for {wildcards.pathogen} cluster {wildcards.cluster}."
    shell:
        "(gatk SelectVariants --reference {input.ref} --variant {input.vcf} "
        "--output {output} --select-type-to-include SNP) 2> {log}"


rule filter_variants:
    input:
        vcf="gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_joint_raw.vcf.gz",
        ref=get_ref_fasta,
    log:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_filtered.log",
    output:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_filtered.vcf.gz",
    message:
        "Filtering the VCF file for {wildcards.pathogen} cluster {wildcards.cluster}."
    shell:
        "(gatk VariantFiltration --reference {input.ref} --variant {input.vcf} --output {output} "
        '--filter-name "Q_and_DP_filter" --filter-expression "QUAL >= 30.0 && DP >= 5") 2> {log}'


rule snp_table:
    input:
        vcf="gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_filtered.vcf.gz",
        ref=get_ref_fasta,
    log:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_filter_table.log",
    output:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_filter_table.tsv",
    message:
        "Outputting the filtered SNPs into a table for {wildcards.pathogen} cluster {wildcards.cluster}."
    shell:
        "(gatk VariantsToTable --variant {input.vcf} --reference {input.ref} --output {output} "
        "--fields CHROM --fields POS --fields ID --fields QUAL --fields REF --fields ALT --fields AC "
        "--fields AF --fields DP --fields HET --fields HOM-REF --fields HOM-VAR --fields NO-CALL --fields VAR "
        "--fields NSAMPLES --fields NCALLED --fields MULTI-ALLELIC --split-multi-allelic "
        "--show-filtered) 2> {log}"


rule sfs:
    input:
        vcf="gatk_{pathogen}/{pathogen}_cluster_{cluster}_SNP_filtered.vcf.gz",
    log:
        "gatk_{pathogen}/{pathogen}_cluster_{cluster}_sfs.log",
    output:
        table="gatk_{pathogen}/{pathogen}_cluster_{cluster}_sfs.tsv",
        barplot="gatk_{pathogen}/{pathogen}_cluster_{cluster}_sfs.pdf",
    message:
        "Calculating the SFS for {wildcards.pathogen} cluster {wildcards.cluster}."
    shell:
        "(Rscript scripts/sfs.R {input} {output.barplot} {output.table}) &> {log}"
