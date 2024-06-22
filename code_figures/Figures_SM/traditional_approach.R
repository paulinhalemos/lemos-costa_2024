## Code showing the plots for the traditional approach, based on regression
 # Steps: 1. Simulate an ultrametric tree with nsp species
 #        2. Get the variance-covariance matrix describing the tree (or a dist matrix)
 #        3. Build a GLV system using info from the tree - same growth rate for everyone
 #        4. Fit the traditional models and see if they recover phylogenetic signal
rm(list = ls())
library(ape)
library(tidyverse)
library(ggpubr)
library(grid)
library(gridExtra)
library(phytools)
library(adephylo)
library(picante)
# 1. Building the ultrametric tree
set.seed(15)
nsp <- 10
spnames <- paste("sp", 1:nsp, sep = "_")
# This function already creates an ultrametric tree (pbtree from phytools does the same)
tree <- rtree(n = nsp, tip.label = spnames)
# 2. Get the variance-covariance matrix of a tree
sigma <- vcv.phylo(tree)

A <- sigma + mean(sigma)

# Calculating Faith's PD (measure of phylogenetic distance)
PD <- sum(tree$edge.length)

# 3. Get all feasible equilibria
results <- data.frame()

zeros <- rep(0, nsp)
ones <- rep(1, nsp)

for (i in 1:(2^nsp - 1)){
  present <- as.numeric(intToBits(i)[1:nsp]) > 0
  k <- sum(present)
  onesk <- ones[1:k]
  composition <- paste(spnames[present], collapse = "-")
  Ak <- A[present, present, drop = FALSE]
  xstar <- solve(Ak, onesk)
  x <- zeros
  x[present] <- xstar
  biomass <- NA
  phylo_comm <- t(as.matrix(present*1))
  colnames(phylo_comm) <- spnames
  if (all(xstar > 0)) {
    tmp <- data.frame(t(x))
    colnames(tmp) <- spnames
    biomass <- sum(xstar)
    phylo_dist <- pd(phylo_comm, tree, include.root = TRUE)
    results <- rbind(results, cbind(tmp, data.frame(community = i, k = k, n = nsp, 
                                                    composition = composition, biomass = biomass,
                                                    phylo = phylo_dist)))
  }
}


# Add noise to the data
noise <- rnorm(length(results$biomass))*0.1
#noise <- 0

y_og <- results$biomass
y <- y_og+noise

# 4. Fit the traditional models
# Model with species diversity + phylogenetic distance
m1 <- lm(biomass ~ k + phylo.PD, data = results)
summary(m1)

# Species richness and PD interact
m2 <- lm(biomass ~ k + phylo.PD + k : phylo.PD, data = results)
summary(m2)

## Including presence absence
P <- data.frame((results[,1:nsp]>0)*1) # design matrix
par(mfrow=c(2,3))
# Just presence/absence
m3 <- lm(y ~ (.), data = P)
(cor_m3 <- round(cor(y, m3$fitted.values), 4))
npars_m3 <- length(coef(m3))
plot(m3$fitted.values, y, main = c("y~P"), ylab = "obs biomass", xlab = "predicted"); abline(0,1)
mtext(paste("corr =", cor_m3, " n_pars =", npars_m3))
# Presence/absence with phylogenetic distance
m4 <- lm(y ~ (.) + results$phylo.PD, data = P)
(cor_m4 <- round(cor(y, m4$fitted.values), 4))
npars_m4 <- length(coef(m4))
plot(m4$fitted.values, y, main = c("y~P+phylo_dist"), ylab = "obs biomass", xlab = "predicted"); abline(0,1)
mtext(paste("corr =", cor_m4, " n_pars =", npars_m4))
# Presence/absence with phylogenetic distance and their interaction
m6 <- lm(y ~ (.) + results$phylo.PD + (.):results$phylo.PD, data = P)
(cor_m6 <- round(cor(y, m6$fitted.values), 4))
npars_m6 <- length(coef(m6))
plot(m6$fitted.values, y, main = c("y~P+phylo+(P:phylo)"), ylab = "obs biomass", xlab = "predicted"); abline(0,1)
mtext(paste("corr =", cor_m6, " n_pars =", npars_m6))
# Presence/absence and pairwise
m3a <- lm(y ~ (.)^2, data = P)
(cor_m3a <- round(cor(y, m3a$fitted.values), 4))
npars_m3a <- length(coef(m3a))
plot(m3a$fitted.values, y, main = c("y~P+P_pair"), ylab = "obs biomas", xlab = "predicted"); abline(0,1)
mtext(paste("corr =", cor_m3a, " n_pars =", npars_m3a))
# Presence/absence with phylogenetic distance and all pairwise
m5 <- lm(y ~ (.)^2 + results$phylo.PD, data = P)
cor_m5 <- round(cor(y, m5$fitted.values), 4)
npars_m5 <- length(coef(m5))
plot(m5$fitted.values, y, main = c("y~P+P_pair+phylo_dist"), ylab = "obs biomass", xlab = "predicted"); abline(0,1)
mtext(paste("corr =", cor_m5, " n_pars =", npars_m5))
# Presence/absence with phylogenetic distance, all pairwise and interaction between presence and dist
m7 <- lm(y ~ (.)^2 + results$phylo.PD + (.):results$phylo.PD, data = P)
(cor_m7 <- round(cor(y, m7$fitted.values), 4))
npars_m7 <- length(coef(m7))
plot(m7$fitted.values, y, main = c("y~P+P_pair+phylo+(P:phylo)"), ylab = "obs biomass", xlab = "predicted"); abline(0,1)
mtext(paste("corr =", cor_m7, " n_pars =", npars_m7))

par(mfrow=c(1,1))


## Function to plot the diagnostic plots using ggplot
diagPlot<-function(model, obs){
  p0<-ggplot(model, aes(.fitted, obs))+geom_point(alpha = 0.5)
  p0<-p0+stat_smooth(method="lm", alpha = 0.5)
  p0<-p0+xlab("Fitted values")+ylab("Observed")
  p0<-p0+ggtitle("Observed vs Fitted Plot")+theme_bw()
  
  p1<-ggplot(model, aes(.fitted, .resid))+geom_point(alpha = 0.5)
  p1<-p1+stat_smooth(method="loess")+geom_hline(yintercept=0, col="red", linetype="dashed")
  p1<-p1+xlab("Fitted values")+ylab("Residuals")
  p1<-p1+ggtitle("Residual vs Fitted Plot")+theme_bw()
  
  p2<-ggplot(model, aes(qqnorm(.stdresid)[[1]], .stdresid))+geom_point(na.rm = TRUE, alpha = 0.5)
  p2<-p2+geom_abline()+xlab("Theoretical Quantiles")+ylab("Standardized Residuals")
  p2<-p2+ggtitle("Normal Q-Q")+theme_bw()
  
  p3<-ggplot(model, aes(.fitted, sqrt(abs(.stdresid))))+geom_point(na.rm=TRUE, alpha = 0.5)
  p3<-p3+stat_smooth(method="loess", na.rm = TRUE)+xlab("Fitted Value")
  p3<-p3+ylab(expression(sqrt("|Standardized residuals|")))
  p3<-p3+ggtitle("Scale-Location")+theme_bw()
  
  p4<-ggplot(model, aes(seq_along(.cooksd), .cooksd))+geom_bar(stat="identity", position="identity")
  p4<-p4+xlab("Obs. Number")+ylab("Cook's distance")
  p4<-p4+ggtitle("Cook's distance")+theme_bw()
  
  p5<-ggplot(model, aes(.hat, .stdresid))+geom_point(aes(size=.cooksd), na.rm=TRUE, alpha = 0.7)
  #p5<-p5+stat_smooth(method="loess", na.rm=TRUE)
  p5<-p5+xlab("Leverage")+ylab("Standardized Residuals")
  p5<-p5+ggtitle("Residual vs Leverage Plot")
  p5<-p5+scale_size_continuous("Cook's Distance", range=c(1,5))
  p5<-p5+theme_bw()+theme(legend.position="bottom")
  
  return(list(predPlot=p0, rvfPlot=p1, qqPlot=p2, sclLocPlot=p3, cdPlot=p4, rvlevPlot=p5))
}

# For model 3
diagPlts_3<-diagPlot(m3, y)
plot_m3 <- ggarrange(diagPlts_3$predPlot, diagPlts_3$rvfPlot, diagPlts_3$qqPlot, diagPlts_3$sclLocPlot, ncol=2, nrow = 2)
plotm0 <- annotate_figure(plot_m3, top = text_grob(paste("M0: y~P ", " corr=", cor_m3, " n_pars=", npars_m3), size = 14))
plotm0
#ggsave("traditional_m0.jpg", plot = plotm0)
# For model 4
diagPlts_4<-diagPlot(m4, y)
plot_m4 <- ggarrange(diagPlts_4$predPlot, diagPlts_4$rvfPlot, diagPlts_4$qqPlot, diagPlts_4$sclLocPlot, ncol=2, nrow = 2)
plotm1 <- annotate_figure(plot_m4, top = text_grob(paste("M1: y~P+phylo_dist ", " corr=", cor_m4, " n_pars=", npars_m4), size = 14))
plotm1
#ggsave("traditional_m1.jpg", plot = plotm1)

# For model 6
diagPlts_6<-diagPlot(m6, y)
plot_m6 <- ggarrange(diagPlts_6$predPlot, diagPlts_6$rvfPlot, diagPlts_6$qqPlot, diagPlts_6$sclLocPlot, ncol=2, nrow = 2)
plotm2 <- annotate_figure(plot_m6, top = text_grob(paste("M2: y~P+phylo+(P:phylo) ", " corr=", cor_m6, " n_pars=", npars_m6), size = 14))
plotm2
#ggsave("traditional_m2.jpg", plot = plotm2)

# For model 5
diagPlts_5<-diagPlot(m5, y)
plot_m5 <- ggarrange(diagPlts_5$predPlot, diagPlts_5$rvfPlot, diagPlts_5$qqPlot, diagPlts_5$sclLocPlot, ncol=2, nrow = 2)
plotm3 <- annotate_figure(plot_m5, top = text_grob(paste("M3: y~P+phylo_dist ", " corr=", cor_m5, " n_pars=", npars_m5), size = 14))
plotm3
#ggsave("traditional_m3.jpg", plot = plotm3)

# For model 7
diagPlts_7<-diagPlot(m7, y)
plot_m7 <- ggarrange(diagPlts_7$predPlot, diagPlts_7$rvfPlot, diagPlts_7$qqPlot, diagPlts_7$sclLocPlot, ncol=2, nrow = 2)
plotm4 <- annotate_figure(plot_m7, top = text_grob(paste("M4: y~P+P_pair+phylo+(P:phylo) ", " corr=", cor_m7, " n_pars=", npars_m7), size = 14))
plotm4
#ggsave("traditional_m4.jpg", plot = plotm4)
