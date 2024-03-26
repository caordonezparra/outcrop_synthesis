library(targets)
library(readr)
library(dplyr)

source("create_barplot.R")
source("figure_panel.R")

tar_option_set(packages = c("readr", "dplyr", "ggplot2", "cowplot", "sv"))

prep_veg <- function(fdat){
  vegetation <- fdat %>%
    mutate(Vegetation = as.factor(Vegetation),
           dummy = "dummy")
  
  ## Now, we'll indicate the order that we want each vegetation type to appear.
  
  vegetation$Vegetation <- ordered(vegetation$Vegetation,
                                   levels = c("Inselberg",
                                              "Campo de altitude",
                                              "Canga",
                                              "Campo rupestre"))
  return(vegetation)
}

vegetationColors <- c("Inselberg" = "#DECC63", "Campo de altitude" = "#799664",
                      "Canga" = "#997054", "Campo rupestre" = "#B3B4AE")

list(
  tar_target(raw_myfile, "/Volumes/Apollo M100/2. Doutorado/Projects/outcrop_synthesis/vegetation.csv", format = "file"),
  tar_target(mydat, read_csv(raw_myfile)),
  tar_target(vegetation, prep_veg(mydat)),
  tar_target(Figure_1E1, create_barplot(vegetation, nStudies, "Studies Percentage (%)")),
  tar_target(Figure_1E2, create_barplot(vegetation, nSpecies, "Taxa Percentage (%)")),
  tar_target(Figure_1E, figure_panel(Figure_1E1, Figure_1E2)),
  tar_target(name = Figure_1E_tofile, command = ggsave("Figure_1E.svg", plot = Figure_1E, width = 17.25, height = 5.5, units = "cm", dpi = 600), format = "file"))