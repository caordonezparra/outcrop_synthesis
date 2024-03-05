library(targets)
library(readr)
library(dplyr)

source("create_barplot.R")

tar_option_set(packages = c("readr", "dplyr", "ggplot2"))

vegetation <- read_csv("/Volumes/Apollo M100/2. Doutorado/Projects/outcrop_synthesis/vegetation.csv")

vegetation <- vegetation %>%
  mutate(Vegetation = as.factor(Vegetation),
         dummy = "dummy")

# Now, we'll indicate the order that we want each vegetation type to appear.

vegetation$Vegetation <- ordered(vegetation$Vegetation,
                                 levels = c("Inselberg",
                                            "Campo de altitude",
                                            "Canga",
                                            "Campo rupestre"))

vegetationColors <- c("Inselberg" = "#DECC63", "Campo de altitude" = "#799664",
                      "Canga" = "#997054", "Campo rupestre" = "#B3B4AE")

list(
  tar_target(Figure_1E1, create_barplot(vegetation, nStudies, "Studies Percentage (%)")),
  tar_target(Figure_1E2, create_barplot(vegetation, nSpecies, "Taxa Percentage (%)"))
)