## Molecular vs fitted branch lengths
 # Code to compare the predicted biomasses from molecular branch lengths
 # And compare the fitteed branch lengths and the molecular branch lengths
rm(list = ls())

library(tidyverse)
library(ggpubr)

source("functions_aux.R")



## Biodiversity II ordering was done manually
# Load the tree
data_sets <- c("bioII_2001-6", "bioII_2006-7", "bioII_2007-8", "bioII_2008-7",
               "bioII_2011-8", "bioII_2012-7", "bioII_2014-8", "bioII_2015-8", "bioII_2017-7")

data_index <- 9

load(paste0("../data/biodiversity_II/organized_data/", data_sets[data_index], "_tree.RData"))
load(paste0("../organized_results/", data_sets[data_index], "_1.RData"))
# Read the empircal data
Obs <- t(as.matrix(read.csv(paste0("../data/biodiversity_II/organized_data/", data_sets[data_index], ".csv"))))
# V matrix associated with the tree
V <- as.matrix(read.csv(paste0("../data/biodiversity_II/organized_data/", data_sets[data_index], "_V.csv")))




# Compact version of data - speeds up calculations
P <- (Obs > 0) * 1
# Number of species in each community
nsp <- colSums(P)
# Parameters
n <- nrow(P)

# Here we had to do the new order manually
all_ordering <- list("Bio2001" = c(1,3,6,9,13,14,15,7,16,10,17,18,2,4,5,19,8,11,20,21,22,23,12,24,25),
                     "Bio2006" = c(1,3,5,9,20,13,21,22, 14,23,24, 7,25,10,15,26,27,28, 2,4,6,16,29,30, 11,17,31,32,33, 18,34,35, 8,12,36,19,37,38,39),
                     "Bio2007" = c(1,3,5,9,20,13,21,22, 14,23,24, 7,25,10,15,26,27,28, 2,4,6,16,29,30, 11,17,31,32,33, 18,34,35, 8,12,36,19,37,38,39),
                     "Bio2008" = c(1,3,6,8,17,11,18,19,20, 7,21,9,12,22,23,24, 2,4,5,13,25,26, 10,14,27,28,29, 15,30,31,16,32,33),
                     "Bio2011" = c(1,2,5,6,9,20,13,21,22, 14,23,24,25, 8,26,10,15,27,28,29, 3,4,7,16,30,31, 11,17,32,33,34, 18,35,36,12,37,19,38,39),
                     "Bio2012" = c(1,2,6,8,16,11,17,18,19, 7,20,9,12,21,22,23, 3,4,5,13,24,25, 10,14,26,27,28,29, 15,30,31),
                     "Bio2014" = c(1,2,6,8,16,11,17,18,19, 7,20,9,12,21,22,23, 3,4,5,13,24,25, 10,14,26,27,28,29, 15,30,31),
                     "Bio2015" = c(1,2,4,6,9,18,12,19,20, 13,21,22,23, 8,24,10,14,25,26,27, 3,5,7,15,28,29, 11,16,30,31,32,33, 17,34,35),
                     "Bio2017" = c(1,2,4,6,9,19,13,20,21, 14,22,23,24, 8,25,10,15,26,27,28, 3,5,7,16,29,30, 11,17,31,32,33,34, 12,35,18,36,37))


new_order <- all_ordering[[data_index]]

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
ggsave(filename = paste0("figs_mol_x_fit/", "molecular_fit_", results$label, ".pdf"))
