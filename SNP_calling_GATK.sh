### Merge the files ###

samtools merge -f merged_bams/OL1008_adRm_sorted_rmdup_m.bam capture-data/OL1008_adRm_sorted_rmdup.bam
samtools merge -f merged_bams/OL1099_adRm_sorted_rmdup_m.bam capture-data/OL1099_adRm_sorted_rmdup.bam
samtools merge -f merged_bams/OL1128_adRm_sorted_rmdup_m.bam capture-data/OL1128_adRm_sorted_rmdup.bam
samtools merge -f merged_bams/OL1214_adRm_sorted_rmdup_m.bam capture-data/OL1214_adRm_sorted_rmdup.bam
samtools merge -f merged_bams/OL1231_adRm_sorted_rmdup_m.bam capture-data/OL1231_adRm_sorted_rmdup.bam
samtools merge -f merged_bams/OL1385_adRm_sorted_rmdup_m.bam capture-data/OL1385_adRm_sorted_rmdup.bam
samtools merge -f merged_bams/OL1389_adRm_sorted_rmdup_m.bam capture-data/OL1389_adRm_sorted_rmdup.bam capture_second_round/OL1389_sorted_rmdup.bam capture_second_round/OL1389ii_sorted_rmdup.bam
samtools merge -f merged_bams/OL1584_adRm_sorted_rmdup_m.bam capture-data/OL1584_adRm_sorted_rmdup.bam capture_second_round/OL1584_sorted_rmdup.bam capture_second_round/OL1584ii_sorted_rmdup.bam
samtools merge -f merged_bams/OL1585_adRm_sorted_rmdup_m.bam capture-data/OL1585_adRm_sorted_rmdup.bam
samtools merge -f merged_bams/OL1599_adRm_sorted_rmdup_m.bam capture-data/OL1599_adRm_sorted_rmdup.bam capture_second_round/OL1599_sorted_rmdup.bam capture_second_round/OL1599ii_sorted_rmdup.bam
samtools merge -f merged_bams/OL1934_adRm_sorted_rmdup_m.bam capture-data/OL1934_adRm_sorted_rmdup.bam capture_second_round/OL1934_sorted_rmdup.bam capture_second_round/OL1934ii_sorted_rmdup.bam
samtools merge -f merged_bams/OL1936_adRm_sorted_rmdup_m.bam capture-data/OL1936_adRm_sorted_rmdup.bam capture_second_round/OL1936_sorted_rmdup.bam capture_second_round/OL1936ii_sorted_rmdup.bam
samtools merge -f merged_bams/OL1984_adRm_sorted_rmdup_m.bam capture-data/OL1984_adRm_sorted_rmdup.bam capture_second_round/OL1984_sorted_rmdup.bam capture_second_round/OL1984ii_sorted_rmdup.bam
samtools merge -f merged_bams/OL1986_adRm_sorted_rmdup_m.bam capture-data/OL1986_adRm_sorted_rmdup.bam capture_second_round/OL1986_sorted_rmdup.bam capture_second_round/OL1986ii_sorted_rmdup.bam
samtools merge -f merged_bams/OL1987_adRm_sorted_rmdup_m.bam capture-data/OL1987_adRm_sorted_rmdup.bam
samtools merge -f merged_bams/OL1999_adRm_sorted_rmdup_m.bam capture-data/OL1999_adRm_sorted_rmdup.bam
samtools merge -f merged_bams/OL2267_adRm_sorted_rmdup_m.bam capture_second_round/OL2267ii_sorted_rmdup.bam capture_second_round/OL2267iii_sorted_rmdup.bam
samtools merge -f merged_bams/OL2268_adRm_sorted_rmdup_m.bam capture_second_round/OL2268i_sorted_rmdup.bam capture_second_round/OL2268ii_sorted_rmdup.bam capture_second_round/OL2268iii_sorted_rmdup.bam
samtools merge -f merged_bams/OL2272_adRm_sorted_rmdup_m.bam capture_second_round/OL2272i_sorted_rmdup.bam capture_second_round/OL2272ii_sorted_rmdup.bam capture_second_round/OL2272iii_sorted_rmdup.bam capture_second_round/OL2272_merged_rmdup_q20.bam
samtools merge -f merged_bams/OL2178_adRm_sorted_rmdup_m.bam capture_second_round/OL2178i_sorted_rmdup.bam capture_second_round/OL2178ii_sorted_rmdup.bam capture_second_round/OL2178iii_sorted_rmdup.bam



### SNP CALLLING NEW ### 


# Trim ends of the reads - 5 bp from each side - from the BAM files with bamUtil

for i in ./merged_bams/*bam; do ~/bin/bamUtil/bin/bam trimBam $i trimmed_bams_snp/`basename -s .bam $i`_trim.bam -L 5 -R 5; done


# Fix the Read Group identifiers in Laurent's BAM files with Picard tools 

for i in trimmed_bams_snp/*.bam; do java -jar ~/bin/picard.jar AddOrReplaceReadGroups I=$i O=trimmed_bams_snp/`basename -s .bam $i`_RG.bam RGID=M_A00363:131 RGLB=lib1 RGPL=Illumina RGPU=HHY2CDSXX:2 RGSM=`basename -s _adRm_sorted_rmdup_m_trim.bam $i`; done 


# index bam files for SNP calling 

for i in trimmed_bams_snp/*RG.bam; do samtools index $i $i.bai; done  

# SNP calling from the BAM files with GATK 


# Joint SNP calling - every sample separate 

# Step two from the gatk best practices: https://software.broadinstitute.org/gatk/documentation/article.php?id=3893 and https://software.broadinstitute.org/gatk/documentation/article?id=11813 and https://software.broadinstitute.org/gatk/documentation/article.php?id=4150

# do it with parallel 

parallel --jobs 5 '~/bin/gatk-4.1.2.0/gatk HaplotypeCaller --input trimmed_bams_snp/{1}.bam --output joint_SNP_calling/{1}.g.vcf --reference ref_genome/EF523390.1_mask.fasta --min-base-quality-score 30 --output-mode EMIT_ALL_SITES --sample-ploidy 1 -ERC GVCF --do-not-run-physical-phasing true' ::: OL1008_adRm_sorted_rmdup_m_trim_RG OL1099_adRm_sorted_rmdup_m_trim_RG OL1128_adRm_sorted_rmdup_m_trim_RG OL1214_adRm_sorted_rmdup_m_trim_RG OL1231_adRm_sorted_rmdup_m_trim_RG OL1385_adRm_sorted_rmdup_m_trim_RG OL1389_adRm_sorted_rmdup_m_trim_RG OL1584_adRm_sorted_rmdup_m_trim_RG OL1585_adRm_sorted_rmdup_m_trim_RG OL1599_adRm_sorted_rmdup_m_trim_RG OL1934_adRm_sorted_rmdup_m_trim_RG OL1936_adRm_sorted_rmdup_m_trim_RG OL1984_adRm_sorted_rmdup_m_trim_RG OL1986_adRm_sorted_rmdup_m_trim_RG OL1987_adRm_sorted_rmdup_m_trim_RG OL1999_adRm_sorted_rmdup_m_trim_RG OL2178_adRm_sorted_rmdup_m_trim_RG OL2267_adRm_sorted_rmdup_m_trim_RG OL2268_adRm_sorted_rmdup_m_trim_RG OL2272_adRm_sorted_rmdup_m_trim_RG



# aggregate all the VCFs into one and then FilterVariants after. Because we've got haploid data we use the CombineGVCFs tool instead of the new GenomicsDBImport

~/bin/gatk-4.1.2.0/gatk CombineGVCFs  --reference ref_genome/EF523390.1_mask.fasta --variant joint_SNP_calling/OL1008_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1099_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1128_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1214_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1231_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1385_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1389_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1584_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1585_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1599_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1934_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1936_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1984_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1986_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1987_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL1999_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL2178_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL2267_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL2268_adRm_sorted_rmdup_m_trim_RG.g.vcf --variant joint_SNP_calling/OL2272_adRm_sorted_rmdup_m_trim_RG.g.vcf --output joint_SNP_calling/ancient_SNPs_cohort.g.vcf


# Feed the aggregated data to GenotypeGVCFs 

~/bin/gatk-4.1.2.0/gatk GenotypeGVCFs  --reference ref_genome/EF523390.1_mask.fasta --variant joint_SNP_calling/ancient_SNPs_cohort.g.vcf --output joint_SNP_calling/ancient_SNPs_joint_gen_raw.vcf


# Filter the joint raw genotypes with the same tool as yesterday

~/bin/gatk-4.1.2.0/gatk VariantFiltration -R ref_genome/EF523390.1_mask.fasta -V joint_SNP_calling/ancient_SNPs_joint_gen_raw.vcf -O joint_SNP_calling/ancient_SNPs_joint_filtered_4x.vcf --filter-expression "QUAL >= 30.0 && DP >= 4" --filter-name "Q_and_DP_filter"

~/bin/gatk-4.1.2.0/gatk VariantFiltration -R ref_genome/EF523390.1_mask.fasta -V joint_SNP_calling/ancient_SNPs_joint_gen_raw.vcf -O joint_SNP_calling/ancient_SNPs_joint_filtered_6x.vcf --filter-expression "QUAL >= 30.0 && DP >= 6" --filter-name "Q_and_DP_filter"


# Produce a summary table for the SNPs

~/bin/gatk-4.1.2.0/gatk VariantsToTable --variant joint_SNP_calling/ancient_SNPs_joint_filtered_4x.vcf -F CHROM -F POS -F ID -F QUAL -F REF -F ALT -F AC -F DP --reference ref_genome/EF523390.1_mask.fasta --output joint_SNP_calling/ancient_SNPs_joint_filtered_4x_table.csv --split-multi-allelic --show-filtered

~/bin/gatk-4.1.2.0/gatk VariantsToTable --variant joint_SNP_calling/ancient_SNPs_joint_filtered_6x.vcf -F CHROM -F POS -F ID -F QUAL -F REF -F ALT -F AC -F DP --reference ref_genome/EF523390.1_mask.fasta --output joint_SNP_calling/ancient_SNPs_joint_filtered_6x_table.csv --split-multi-allelic --show-filtered



# I trust the joint SNP calling as it gives me similar results to my bcftools SNP calling pipeline 


### PILE UP SEQUENCES ###

# We use the HTSBOX PILEUP Tool - trim 5 bp from each end of the reads and drop any allele that has a depth of less than 4 reads 


# From the non trim bam files 


for i in merged_bams/*.bam; do echo $i; /home/alex/bin/htsbox pileup -f ref_genome/EF523390.1_mask.fasta -l15 -T5 -q30 -Q30 -M -s4 $i > pileup_fastas/`basename $i`.fasta; done 

for i in pileup_fastas/*fasta; do sed -i "s/>EF523390.1/>`basename -s _adRm_sorted_rmdup_m.bam.fasta $i`/g" $i; done

for i in pileup_fastas/*fasta; do echo "mv $i pileup_fastas/`basename -s _adRm_sorted_rmdup_m.bam.fasta $i`.fasta"; done 


# Average depth for each pile up fasta seq OLD 


merged_bams/OL1008_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': -nan
merged_bams/OL1099_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 6.67
merged_bams/OL1128_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': -nan
merged_bams/OL1214_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': -nan
merged_bams/OL1231_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.00
merged_bams/OL1385_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 38.22
merged_bams/OL1389_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 7.94
merged_bams/OL1584_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 8.49
merged_bams/OL1585_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.75
merged_bams/OL1599_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 5.61
merged_bams/OL1934_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.89
merged_bams/OL1936_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 5.96
merged_bams/OL1984_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.89
merged_bams/OL1986_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 6.79
merged_bams/OL1987_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.87
merged_bams/OL1999_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.00
merged_bams/OL2178_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 5.56
merged_bams/OL2267_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.49
merged_bams/OL2268_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 5.28
merged_bams/OL2272_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 7.83


# NEW AVERAGE DEPTH

[M::write_fa] average depth for contig 'EF523390.1': -nan
merged_bams/OL1099_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 6.67
merged_bams/OL1128_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': -nan
merged_bams/OL1214_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': -nan
merged_bams/OL1231_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.00
merged_bams/OL1385_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 38.22
merged_bams/OL1389_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 7.94
merged_bams/OL1584_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 8.49
merged_bams/OL1585_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.75
merged_bams/OL1599_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 5.61
merged_bams/OL1934_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.89
merged_bams/OL1936_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 5.96
merged_bams/OL1984_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.89
merged_bams/OL1986_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 6.79
merged_bams/OL1987_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.87
merged_bams/OL1999_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.00
merged_bams/OL2178_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 5.56
merged_bams/OL2267_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 4.49
merged_bams/OL2268_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 5.28
merged_bams/OL2272_adRm_sorted_rmdup_m.bam
[M::write_fa] average depth for contig 'EF523390.1': 7.83

