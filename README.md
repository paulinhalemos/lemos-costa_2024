# Data and code accompanying "Phylogeny structures species' interactions in experimental ecological communities" by Lemos-Costa, Miller and Allesina 2024

All scripts use relative paths and assume the working directory is the one where the code is located.

Contents:

- **code** code to fit the statistical models as well as the randomizations; code to organize raw results (removing duplicate randomized trees)
    - get_random_V.R code with auxiliary functions to run the model; randomization the structure of the tree
    - maximization.R code with auxiliary functions of the maximization algorithms (hill climbing, simulated annealing, likelihood)
    - organize_results.R code to organize the results
    - pipeline.R function to run the model and save the results
    - README.Rmd markdown showing an example on how to run the model with the data set from Cadotte 

- **data** raw data from original publications; code to parse and organize the raw data
    - README.md explains the structure of the folders
    - biodiversity_II 
    - cadotte_2013
    - van_ruijven
    - general_code folder with a code (process_dataset.R) common to processing all data sets
    

- **results** resulting model fits and randomizations; typically, each model/data set is fitted several times, including about 100 randomizations (will take about 1 day to run on a recent computer); by running several times in parallel, with different random seeds, one can collect a large number of randomizations. As an example, many files containing randomizations for Cadotte (2013) model 1 are provided.

- **functions** stan functions used to speed up code

- **organized_results** organized results for each data set/model combination, obtained by merging several files containing randomizations and removing duplicate randomized trees.

- **code_figures** code to reproduce the figures and tables in the manuscript.
  - table_results.R and build_tibbles_for_figures.R builds tables and tibbles with the results in the appropriate format for plotting
  - lookup_names.csv file with the label associated with each data set (first column) and their short name for plotting (second column)
  - table_results.csv summarizes the results for plotting
  - forplots.RData summarizes the likelihoods and observations in the appropriate format for plot
  - Figure2.R plots figure 2 of the main manuscript 
  - Figure3.R plots figure 3 of the main manuscript
  - Figure4.R plots figure 4 of the main manuscript
  - Folder Figures_SM with code to reproduce the figures in the Supplementary Material



