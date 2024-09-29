# This is the _targets.R file with the pipeline that is going to allow us to create Figure S2.

# The following are the packages that we will need to use to run the whole pipeline

library(targets) # v. 1.5.1 for creating and running pipelines
library(gert) # v. 2.0.1 for managing Git repositories
library(readr) # v. 2.1.5 for loading datasets.
library(dplyr) # v. 1.1.4 for data frame manipulation.
library(tidyr) # v. 1.3.1 for for more database manipulation.
library(plotrix) # v. 3.8-4 for standard error calculation.
library(ggplot2) # v. 3.4.4 for creating different visualizations
library(scales) # v. 1.3.0 for percent scale.
library(cowplot) # v. 1.1.3 for creating plot grids.
library(svglite) # v. 2.1.3 for creating '.svg' files.

# We'll source the functions that we have written to create this figure.

source("joined_datasets.R")
source("process_and_plot.R")
source("figureS2_panel1.R")
source("figure_dormancy.R")
source("figure_syndrome.R")
source("figure_season.R")
source("figureS2_panel2.R")
source("figureS2_panel_final.R")


# To start the pipeline, we will call explicitly the packages the pipeline is going to use.

tar_option_set(packages = c("gert", "stringr", "readr", "dplyr", "tidyr",
                            "ggplot2", "cowplot", "svglite", "plotrix"))

selected_families <- c("Asteraceae", "Eriocaulaceae", "Xyridaceae", "Poaceae", "Cyperaceae", "Bromeliaceae",
                       "Velloziaceae", "Verbenaceae", "Fabaceae", "Melastomataceae")

scale_color_families <- c("All" = "grey", 
                          "Melastomataceae" = "#F9E856", 
                          "Fabaceae" = "#C4DE50",
                          "Velloziaceae" = "#8FCF63", 
                          "Asteraceae" = "#69BD78", 
                          "Bromeliaceae" = "#52A686",
                          "Eriocaulaceae" = "#478F8B", 
                          "Verbenaceae" = "#41768B", 
                          "Xyridaceae" = "#3F5E89",
                          "Poaceae" = "#424484", 
                          "Cyperaceae" = "#432671")

color_dormancy <- c("ND" = "#EF7E87", "D" = "#A95DA2", "Both" = "#D26691")

color_syndrome <- c("Anemochory" = "#FDC300", "Autochory" = "#6FAE58", "Zoochory" = "#EA4743")

color_season <- c("ED" = "#FB8300", "ER" = "#28A6A5", "LD" = "#FFA849", "LR" = "#53D5D4", 
                  "More than 1" = "#F2F2F2")

# Now, we're setting the targets

list(
  tar_target(raw_myfile,  str_c(git_find(),"/Preprocessing/traits_joined.csv"), format = "file"),
  tar_target(mydat, read_csv(raw_myfile)),
  tar_target(traits_all, joined_datasets(mydat)),
  tar_target(FigureS2A_1, process_and_plot(traits_all, "Dry_mass")),
  tar_target(FigureS2A_2, process_and_plot(traits_all, "Water_content")),
  tar_target(FigureS2A_3, process_and_plot(traits_all, "Embryoless")),
  tar_target(FigureS2A_4, process_and_plot(traits_all, "Viable")),
  tar_target(FigureS2B_1, figure_dormancy(traits_all)),
  tar_target(FigureS2B_2, figure_syndrome(traits_all)),
  tar_target(FigureS2B_3, figure_season(traits_all)),
  tar_target(FigureS2A, figureS2_panel1(FigureS2A_1, FigureS2A_2,FigureS2A_3, FigureS2A_4)),
  tar_target(FigureS2B, figureS2_panel2(FigureS2B_1, FigureS2B_2, FigureS2B_3)),
  tar_target(FigureS2, figureS2_panel_final(FigureS2A, FigureS2B)),
  tar_target(name = Figure_S2_tofile, 
             command = ggsave("FigureS2.svg", plot = FigureS2, width = 16, 
                              height = 22.7, units = "cm", dpi = 600), format = "file"))



