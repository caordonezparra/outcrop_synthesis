# The following code is designed to create Figure S2 of Ordóñez-Parra et al, which is a series of boxplots showing the data distribution of the seven evaluated functional traits and the three families with the larges and lowest trait values (among the 10 most represent families in our database).

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: August 1st, 2024.

#### 0. Load the packages (assuming they are already installed) ####

library(gert) # v. 2.0.1 for managing Git repositories
library(stringr) # v. 1.5.1 for string operation.
library(targets) # v. 1.1.1 for pipelines

#### 1. Defining the working directory ####

setwd(str_c(git_find(),"/Figures/FigureS2"))

#### 2. Making the forest plot ####

tar_make()
