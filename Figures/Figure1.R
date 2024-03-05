# This code was designed to produce the barplots that are part of Figure 1 of Ordóñez-Parra et al.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: March 4th, 2024.

#### 0. Set the working directory ####

setwd("/Volumes/Apollo M100/2. Doutorado/Projects/outcrop_synthesis/Figures")

#### 1. Load the packages (assuming they are already installed) ####

library(targets) # v. 1.5.1 for creating pipelines
library(cowplot) # v. 1.1.2 for arranging plot in grids
library(svglite) # v. 2.1.3 for exporting figures as .svg

#### 3. Making the barplots ####

# Figure 1 aims to show the different vegetation types included in the study as well as the number of publications and species studied on each of them. Figures A-D correspond to photographs that were arranged manually in Adobe Illustrator. The following code was used to elaborate Figure 1E, which is barplots showing the proportion of studies and species from each vegetation type (as a percentage of the total of species and studies). 

# Since we are doing two almost identical barplots, we will define a function (see create_barplot.R) and implement an 'targets' pipeline (see the 'targets.R' file in the folder). This pipeline loads the 'vegetation.csv' dataset, which contains the information on studies and species evaluated on each vegetation types, does some pre-processing and produces a bar plot for each. To implement this pipeline, we only need to run the following line:

tar_make()

# Now, that the pipeline has runned, we will arrange these two plots in a grid.

Figure_1E <- plot_grid(tar_read(Figure_1E1), tar_read(Figure_1E2), # The plots we want to put in the grid.
                       ncol = 1, nrow = 2, # The number of columns and rows we want in our grid, respectively.
                       labels = c("E", NULL), # The labels of our plots. 
                       rel_heights = c(1,1.75)) # The relative height of each plot. Since Figure 1F has the legend, we need to adjust its height so that bars have the same height.

Figure_1E

# Now that we're satisfied with the plot, we're going to save this result as a .svg to add it to the photogrids.

ggsave("Figure_1E.svg", width = 17.25, height = 5.5, units = "cm", dpi = 600)
