# Author script for running every bit of code in the proper order
# For MS: XXXX
# Written by: G.Srednick - 3 May 2022
# Revised: 9 October 2023

# System information:
# This code was run on at 16GB Apple M2 Pro with 12 CPUs


# Instructions. 
# These should be run in order. Each script places elements in the environment that allows for successive functions to be run. 
# The only one that can be skipped (due to exceedingly long run-time) is "5_species_interactions.R". 
# The output from this script ("5_species_interactions.R") is saved in the directory and imported in successive scripts.




# (1) Data curation ==== takes about 2 mins to run
# Data for this step should be independently downloaded from DataOne at https://search.dataone.org/view/doi%3A10.6085%2FAA%2FPISCO_kelpforest.1.6 
source("./Manuscript_scripts/1_data_curation.R")


# (2) Get sites and adjusted over-water distances  ==== takes about 15 mins to run
source("./Manuscript_scripts/2_distance_method.R")


# (3) Get temperature variation for each site ==== takes about 15 mins to run
source("./Manuscript_scripts/3_geo_drivers_synchrony.R")

# (4) Run synchrony calculations and get MDS ===== takes about 20 mins to run.
# Ignore the warnings for metaMDS -- stress is acceptable for visualization
# Warnings ok for betadisper
source("./Manuscript_scripts/4_synchrony_stability.R")

# (5) Run correlations between consumers and resources ==== 
#***** WARNING: DONT RUN ME UNLESS YOU WANT TO WAIT FOR >72 HOURS! *******
# exports from this script are supplied and subsequent analyses (e.g., SEM) can be performed without running this script.
source("./Manuscript_scripts/5_species_interactions")

# (6) Run habitat similarity script
source("./Manuscript_scripts/6_hab_similarity.R")

# (7) Run SEM with all the info from above
source("./Manuscript_scripts/7_sem_network.R")
# Ignore the 3 warnings from SEM regarding variances -- model performance has been verified and confirmed.

# (8) Run main figure generator 
source("./Manuscript_scripts/8_figures.R")

# (9) Build structural equation model figures 
source("./Manuscript_scripts/9_SEM_path_plots.R")

# all figures are saved in "./Manuscript_scripts/MS_figures"

## END ## 