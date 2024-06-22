## Molecular vs fitted branch lengths
 # Code to compare the predicted biomasses from molecular branch lengths
 # And compare the fitted branch lengths and the molecular branch lengths
rm(list = ls())

library(tidyverse)
library(ggpubr)

source("functions_aux.R")

# Pick the data set -- Cadotte or Wageningen 
 # FOr Biodiveristy II there is a separate script
# Cadotte
# Load the tree
load("../../data/cadotte_2013/organized_data/cadotte_2013_tree.RData")
load("../../organized_results/cadotte_2013_1.RData")
# Read the empircal data
Obs <- t(as.matrix(read.csv("../../data/cadotte_2013/organized_data/cadotte_2013.csv")))
# V matrix associated with the tree
V <- as.matrix(read.csv(paste0("../../data/cadotte_2013/organized_data/cadotte_2013_V.csv")))

## Wageningen - 
# Load the tree
load("../data/van_ruijven2/organized_data/vr_1_tree.RData")
load("../organized_results/old_organized/vr_1_1.RData")
# Read the empircal data
Obs <- t(as.matrix(read.csv("../data/van_ruijven2/organized_data/vr_1.csv")))
# V matrix associated with the tree
V <- as.matrix(read.csv(paste0("../data/van_ruijven2/organized_data/vr_1_V.csv")))


# Compact version of data - speeds up calculations
P <- (Obs > 0) * 1
# Number of species in each community
nsp <- colSums(P)
# Parameters
n <- nrow(P)


# First we need to get an expanded version of the tree with the branching patterns
tree_str <- expand_tree(tree)
# Now we need to find the right ordering of the tree comparing with our V
new_order <- get_col_redorder(V, tree_str)

# Reorder the columns of V
new_V <- V[new_order,,drop=FALSE]
# Making sure we have the right structure
tmp <- t(new_V[-1,]) %*% diag(tree$edge.length) %*% new_V[-1,]
tmp2 <- vcv(tree)
# This plot should be the 1:1 line
plot(tmp, tmp2);abline(0,1)



# Building the fitted matrix of interaction
A_res <- t(V) %*% diag(results$pars_orig[1:(2*n-1)]) %*% V
# Getting the predictions from the fitted matrix
predictions <- get_endpoints_r(P, A_res, nsp)
# Getting the predictions from the molecular tree
predictions_molecular <- get_endpoints_r(P, vcv(tree)+results$pars_orig[1], nsp)

## Plotting
# Observed and predicted biomasses
pl_ind <- plot_biomass(results, predictions_molecular)

## Molecular branch lengths against fitted branch length
# Parameters need to be reordered too to make the comparison
fitted_pars <- results$pars_orig[1:(2*n-1)][new_order]
dt_pars <- tibble("fitted branch lengths" = fitted_pars[-1], "molecular branch lengths (sqrt)" = sqrt(tree$edge.length))
label_pars <- paste("branch length comparison", results$label)
pl_pars <- ggplot(dt_pars, aes(x = `molecular branch lengths (sqrt)`, y = `fitted branch lengths`)) +
  geom_point() + theme_bw() + ggtitle(label_pars) + scale_x_log10()
show(pl_pars)


ggarrange(pl_ind, pl_pars, nrow = 2)
#ggsave(filename = paste0("figs_mol_x_fit/", "molecular_fit_", results$label, ".pdf"))
