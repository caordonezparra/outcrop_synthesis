# This is the _targets.R file with the pipeline that is going to allow us to create Figure 5.

# The following are the five packages that we will need to use to run the whole pipeline

library(targets) # v. 1.5.1 for creating and running pipelines
library(gert) # v. 2.0.1 for managing Git repositories
library(readr) # v. 2.1.5 for loading datasets.
library(dplyr) # v. 1.1.4 for data frame manipulation.
library(ggplot2) # v. 3.4.4 for creating different visualizations
library(cowplot) # v. 1.1.3 for creating plot grids.
library(svglite) # v. 2.1.3 for creating '.svg' files.

# We'll source the two functions that we have written to create this figure.

source("prep_temperature.R")
source("forest_plot_temperature.R")
source("figure_panel_temperature.R")
source("prep_alternate.R")
source("forest_plot_alternate.R")
source("figure_panel_alternate.R")
source("figure_panel_final.R")

# To start the pipeline, we will call explicitly the packages the pipeline is going to use.

tar_option_set(packages = c("gert", "stringr", "readr", "dplyr", "ggplot2", "cowplot", "svglite"))

scale_color_temperature <- c("Yes" = "#f8d452", "No" = "#7F7F7F", "NA" = "#7F7F7F")

# Now, we're setting the targets

list(
  tar_target(raw_myfile, 
             str_c(git_find(),"/Figures/Figure5/figure_meta_constant.csv"), 
             format = "file"),
  tar_target(mydat, read_csv(raw_myfile)),
  tar_target(constant_temperature, prep_temperature(mydat)),
  tar_target(Figure5a_10, forest_plot_temperature(constant_temperature, 10)),
  tar_target(Figure5a_15, forest_plot_temperature(constant_temperature, 15)),
  tar_target(Figure5a_20, forest_plot_temperature(constant_temperature, 20)),
  tar_target(Figure5a_30, forest_plot_temperature(constant_temperature, 30)),
  tar_target(Figure5a_35, forest_plot_temperature(constant_temperature, 35)),
  tar_target(Figure5a_40, forest_plot_temperature(constant_temperature, 40)),
  tar_target(Figure_5a, figure_panel_temperature(Figure5a_10, Figure5a_15, Figure5a_20, 
                                                 Figure5a_30, Figure5a_35, Figure5a_40)),
  tar_target(raw_myfile2,
             str_c(git_find(),"/Figures/Figure5/figure_meta_alternate.csv"),
             format = "file"),
  tar_target(mydat2, read_csv(raw_myfile2)),
  tar_target(alternate_temperature, prep_alternate(mydat2)),
  tar_target(Figure5b_25, forest_plot_alternate(alternate_temperature, "25/15")),
  tar_target(Figure5b_30, forest_plot_alternate(alternate_temperature, "30/20")),
  tar_target(Figure_5b, figure_panel_alternate(Figure5b_25, Figure5b_30)),
  tar_target(Figure_5, figure_panel_final(Figure_5a, Figure_5b)),
  tar_target(name = Figure_5_tofile, 
             command = ggsave("Figure_5.svg", plot = Figure_5, width = 18, 
                              height = 22.7, units = "cm", dpi = 600), format = "file"))
