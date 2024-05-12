# The following code is designed to create Figure 4 of Ordóñez-Parra et al, which is a is a forest plot showing the effect of light availability (Light vs. Dark) on the germination percentage  of different groups. Moreover, it also shows the interactive effect of seed mass on such effects.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: May 3rd, 2024.

#### 0. Load the packages (assuming they are already installed) ####

library(stringr)
library(gert) # v. 2.0.1 for managing Git repositories
library(readr) # v. 2.1.5 for loading datasets.
library(dplyr) # v. 1.1.4 for data frame manipulation.
library(ggplot2) # v. 3.4.4 for creating different visualizations
library(cowplot) # v. 1.1.3 for creating plot grids.
library(svglite) # v. 2.1.3 for creating '.svg' files.

#### 1. Defining the working directory ####

setwd(str_c(git_find(),"/Figures/Figure4"))

#### 2. Loading and tidying the dataset ####

# For this figure, we'll be using the dataset "figure_meta_light.csv", which contains information on the effect size of light for the different ecological groups tested.

figure_meta_light <- read.csv("figure_meta_light.csv")

# Ensuring that the moderators appear in the order we want them to.

figure_meta_light$Mods <- factor(figure_meta_light$Mods, 
                                 levels = c("Seed mass",
                                            "Xeric",
                                            "Mesic/Xeric",
                                            "Mesic",
                                            "Microhabitat",
                                            "Widespread",
                                            "Restricted",
                                            "Distribution",
                                            "Shrub",
                                            "Herb",
                                            "Growth form",
                                            "All"))

#### 3. Doing the figure and exporting it ####

figure4 <- ggplot(data = figure_meta_light, aes(y = Mods,
                                                 x = post_mean,
                                                 xmin = Lower,
                                                 xmax = Upper,
                                                 col = Significant)) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.8) +
  geom_point(size = 2.5, pch = 15) +
  geom_errorbarh(height=.1) +
  theme_bw(base_size = 12) +
  scale_x_continuous(limits = c(-2.5, 6)) +
  scale_colour_manual(values = c("#7F7F7F", "#FCC38C")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), legend.position = "none", 
        axis.title.y = element_blank(), axis.title.x = element_text(face = "bold"),
        axis.text = element_text(color = "black")) +
  labs(title = "Light", x = "Effect size")

ggsave(filename = "Figure4.svg", device = "svg", height = 8.4, width = 8.4, 
       units = "cm", dpi = 600)
