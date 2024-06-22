# Code to generate Figure 2
source("build_tibbles_for_figures.R")


selected_datasets <- factor(c("Cadotte", "Bio II 7/08", "Wag 03"), levels = c("Cadotte", "Bio II 7/08", "Wag 03"))
pl1 <- plot_histograms_facet_label(likelihoods %>% filter(model == 1, label %in% selected_datasets))

pl2 <- plot_predictions_facet_label(observations %>% filter(model == 1, label %in% selected_datasets))

show(ggarrange(pl2, pl1,labels = c("a)", "b)")))
