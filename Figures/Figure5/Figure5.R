# The following code is designed to create Figure 5 of Ordóñez-Parra et al, which is a is a forest plot showing the effect of different temperatures (both constant and alternate) on the germination percentage and media germination time of different groups. Moreover, it also shows the interactive effect of seed mass on such effects.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: May 3rd, 2024.

#### 0. Load the packages (assuming they are already installed) ####

library(gert) # v. 2.0.1 for managing Git repositories
library(stringr)
library(targets)

#### 1. Defining the working directory ####

setwd(str_c(git_find(),"/Figures/Figure5"))

#### 2. Making the forest plot ####

tar_make()
