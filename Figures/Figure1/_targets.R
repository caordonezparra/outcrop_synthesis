# This is the _targets.R file with the pipeline that is going to allow us to create Figure 1.

# The following are the five packages that we will need to use to run the whole pipeline

library(targets) # v. 1.5.1 for creating and running pipelines
library(gert) # v. 2.0.1 for managing Git repositories
library(readr) # v. 2.1.5 for loading datasets.
library(dplyr) # v. 1.1.4 for data frame manipulation.
library(ggplot2) # v. 3.4.4 for creating different visualizations
library(cowplot) # v. 1.1.3 for creating plot grids.
library(svglite) # v. 2.1.3 for creating '.svg' files.

# We'll source the two functions that we have written to create this figure.

source("prep_veg.R")
source("create_barplot.R") # The one that will create the barplots.
source("figure_panel.R") # The one that will put them together in a grid.

# To start the pipeline, we will call explicitly the packages the pipeline is going to use.

tar_option_set(packages = c("gert", "stringr", "readr", "dplyr", "ggplot2", "cowplot", "svglite"))


vegetationColors <- c("Inselberg" = "#DECC63", "Campo de altitude" = "#799664",
                      "Canga" = "#997054", "Campo rupestre" = "#B3B4AE")

list(
  tar_target(raw_myfile, str_c(git_find(),"/vegetation.csv"), format = "file"),
  tar_target(mydat, read_csv(raw_myfile)),
  tar_target(vegetation, prep_veg(mydat)),
  tar_target(Figure_1E1, create_barplot(vegetation, nStudies, "Studies Percentage (%)")),
  tar_target(Figure_1E2, create_barplot(vegetation, nSpecies, "Taxa Percentage (%)")),
  tar_target(Figure_1E, figure_panel(Figure_1E1, Figure_1E2)),
  tar_target(name = Figure_1E_tofile, command = ggsave("Figure_1E.svg", plot = Figure_1E, width = 17.25, height = 5.5, units = "cm", dpi = 600), format = "file"))
