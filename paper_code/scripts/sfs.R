# Title     : run_fastbaps
# Objective : Cluster genomes with fastbaps
# Created by: Evangelos A. Dimopoulos"
# Created on: 01/12/2021

.libPaths(R.home("/media/jbod/home/antony/bin/miniconda2/envs/mdv/lib/R/library"))

# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
variant_file <- args[1]
sfs_plot <- args[2]
sfs_table <- args[3]

#load libs
library(pegas)
library(ggplot2)


sfs <- function(snps, barplot, out_file) {
  # read the vcf
  vcf <- read.vcf(snps)

  # calculate the sfs
  sfs <- site.spectrum(vcf, folded=TRUE, col='blue')

  sites <- as.vector(sfs)
  samples <- 1:length(sites)
  spectrum <- data.frame(sites, samples)

  # plot the sfs

  sfs_barplot <- ggplot(data=spectrum, aes(x=samples, y=sites)) +
    geom_bar(stat="identity", fill="deepskyblue") +
    theme_minimal() + xlab("SNP frequency") + ylab("Number of SNPs")

  pdf(barplot, paper="a4")
  print(sfs_barplot)
  dev.off()

  # write the spectrum
  write.table(spectrum, file=out_file, row.names=FALSE, quote=FALSE, sep='\t')
}

# run the clustering function
sfs(variant_file, sfs_plot, sfs_table)
