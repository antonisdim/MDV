setwd("~/Marek/second_round_captures")

library(ape)

# outgroup NC_002641.1

mle_best_tree_rooted <- read.nexus('BEAST_root_MLE.nexus')

plot(mle_best_tree_rooted)

plot(drop.tip(mle_best_tree_rooted, c("'NC_002641.1'"), rooted = TRUE, collapse.singles = TRUE))

mle_best_tree_no_out <- drop.tip(mle_best_tree_rooted, c("'NC_002641.1'"), 
                                 rooted = FALSE, collapse.singles = TRUE)

write.nexus(mle_best_tree_no_out, file='BEAST_root_MLE_no_out.nexus') 