For each of three data sets, the corresponding folder contains
- raw data
- code to parse and organize data - scripts named parse_DATASET.R and process_DATASET.R
- organized data used in the manuscript

In all data folders there is a organized_data folder: this folder has three types of files
    1. files that end with *_tree.pdf are plots of the molecular tree
    2. files that end with *_tree.RData refer to the R object with the tree (class phylo)
    3. files that end with *_V.csv refers to the V matrix defining the topology of the tree (see manuscript for details)
    4. files with the data set name with just .csv extension refer to the organized data where columns are species and rows are experiments and the value refers to the biomass of a given species (column) in a given experiment (row)

The folder general_code contains code (process_dataset.R) common to all three data set pipelines


**biodiversity_II** data from the Biodiversity II experiment 
  - parse_biodiversityII.R code to parse the raw data
  - process_biodiversityII.R organizes the data, builds the tree and saves it in the "organized_data" folder
  - species_list.csv list of species that appear in the raw data and their classification as non_target species or target species
  - Spp_codes.csv has the taxonomic information of the selected species with their name in the column species, their genus, family and code is their normalized names


**cadotte_2013** data from the by Cadotte 2013
  - process_cadotte_2013.R organizes the data, builds the tree and saves it in the "organized_data" folder
  - Spp_codes.csv has the taxonomic information of the selected species with their name in the column species, their genus, family and code is their normalized names


**van_ruijven** data from the Wegeningen experiment
  - process_van_ruijven.R organizes the data, builds the tree and saves it in the "organized_data" folder
  - vanruijven-specbiomass-2000-2008.xlsx contains the raw data
  - species_codes.csv has the taxonomic information of the selected species with their name in the column species, their genus, family and code is their normalized names

