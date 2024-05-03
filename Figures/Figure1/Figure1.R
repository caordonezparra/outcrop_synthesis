# This code was designed to produce the barplots that are part of Figure 1 of Ordóñez-Parra et al.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: April 2nd, 2024.

#### 0. Load the packages (assuming they are already installed) ####

library(gert) # v. 2.0.1 for managing Git repositories
library(stringr) # v. 1.5.1 for manipulating strings. 
library(targets) # v. 1.5.1 for creating and running pipelines

#### 1. Set the working directory ####

setwd(str_c(git_find(),"/Figures/Figure1"))

#### 2. Making the barplots ####

# Figure 1 aims to show the different vegetation types included in the study as well as the number of publications and species studied on each of them. Figures A-D correspond to photographs that were arranged manually in Adobe Illustrator. The following code was used to elaborate Figure 1E, which is bar plots showing the proportion of studies and species from each vegetation type (as a percentage of the total of species and studies). 

# Since we are doing two almost identical barplots, we will define a function (see create_barplot.R) and implement an 'targets' pipeline (see the 'targets.R' file in the folder). This pipeline loads the 'vegetation.csv' dataset, which contains the information on studies and species evaluated on each vegetation types, does some pre-processing and produces a bar plot for each. 

# The target pipeline also creates a panel that joins these two barplots in the same figure (see figure_panel.R) and saves it as a .svg file that can be latter manually edited. To implement the whole pipeline, we only need to run the following line:

tar_make()
