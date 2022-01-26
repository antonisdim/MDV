# Title     : run_fastbaps
# Objective : Cluster genomes with fastbaps
# Created by: Evangelos A. Dimopoulos"
# Created on: 01/12/2021

.libPaths(R.home("/media/jbod/home/antony/bin/miniconda2/envs/mdv/lib/R/library"))

# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
input_tree <- args[1]
output_tree <- args[2]

#load libs
library(ape)

prune_tree <- function(input_tree, output_tree) {
  # read the tree
  mle_tree_rooted <- read.tree(input_tree)

  # prune the HVT root
  mle_tree_rooted_no_out <- drop.tip(mle_tree_rooted, c("NC_002641.1"),
                                   rooted = TRUE, collapse.singles = TRUE)

  mle_tree_rooted_no_out_bifur <- multi2di(mle_tree_rooted_no_out, random=TRUE)

  # write the resulting unrooted tree
  write.nexus(mle_tree_rooted_no_out_bifur, file=output_tree)
}

# run the clustering function
prune_tree(input_tree, output_tree)
