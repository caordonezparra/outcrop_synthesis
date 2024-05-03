# The following code is designed to create Figure 6 of Ordóñez-Parra et al, which is a is a forest plot showing the effect of fire-related cues on germination percentage and time of different groups. Moreover, it also shows the interactive effect of seed mass on such effects.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: May 3rd, 2024.

#### 0. Load the packages (assuming they are already installed) ####

library(gert) # v. 2.0.1 for managing Git repositories
library(stringr)
library(readr) # v. 2.1.5 for loading datasets.
library(dplyr) # v. 1.1.4 for data frame manipulation.
library(ggplot2) # v. 3.4.4 for creating different visualizations
library(cowplot) # v. 1.1.3 for creating plot grids.
library(svglite) # v. 2.1.3 for creating '.svg' files.

#### 1. Defining the working directory ####

setwd(str_c(git_find(),"/Figures/Figure6"))

#### 2. Load the database ####

# We are loading the dataset with the information of the effect sizes of the different fire-related cues

figure_meta_fire <-read.csv("figure_meta_fire.csv")

scale_colors_fire <- c("Yes" = "#FB9F83", "No" = "#7F7F7F", "NA" = "#7F7F7F")

#### 3. Elaborate the different figures and join them in a panel ####

# First, the forest plot with all fire-related cues. We will start by filtering the database to include the models that evaluated all cues together.

figure_meta_fire_all <- figure_meta_fire %>%
  filter(Cue == "All")

# Now, we will we tidy up the dataset so that moderators appear in the order we want.

figure_meta_fire_all$Mods <- factor(figure_meta_fire_all$Mods, 
                                    levels = c("Seed mass",
                                               "Non-dormant",
                                               "Dormant",
                                               "Dormancy",
                                               "Smoke",
                                               "Heat",
                                               "Single",
                                               "Combined",
                                               "Treatment type",
                                               "All"))

figure_meta_fire_all$Variable <- factor(figure_meta_fire_all$Variable, 
                                        levels = c("Percentage",
                                                   "Time"))

# With these dataset, we are going to to produce our figure.

figure6a <- ggplot(data = figure_meta_fire_all, aes(y = Mods,
                                                    x = post_mean,
                                                    xmin = Lower,
                                                    xmax = Upper,
                                                    col = Significant)) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.8) +
  geom_point(size = 1.8, shape = 15, position = position_dodge(width = .6)) +
  geom_errorbarh(height= 0, position = position_dodge(width = .6)) +
  theme_bw(base_size = 8) +
  scale_colour_manual(values = scale_colors_fire) +
  xlim(-6, 4) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(), legend.position = "none", 
        axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold")) +
  labs(title = "Fire-related cues", x = "Effect size") +
  facet_wrap(~Variable)

figure6a

# We'll be repeating the same procedure with the analysis of the effect of the different heat shock treatments.

# Filtering the dataset.

figure_meta_fire_shock <- figure_meta_fire %>%
  filter(Cue == "Shock")

# Ensuring the models appear in the order we want.

figure_meta_fire_shock$Mods <- factor(figure_meta_fire_shock$Mods, 
                                           levels = c("200 °C for 1'",
                                                      "100 °C for 5'",
                                                      "100 °C for 1'",
                                                      "80 °C for 5'",
                                                      "80 °C for 2'"))



figure_meta_fire_shock$Variable <- factor(figure_meta_fire_shock$Variable, 
                                          levels = c("Percentage",
                                                     "Time"))

# Creating the figure for the heat shock data.

figure6b <- ggplot(data = figure_meta_fire_shock, aes(y = Mods,
                                                      x = post_mean,
                                                      xmin = Lower,
                                                      xmax = Upper,
                                                      col = Significant)) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.8) +
  geom_point(size = 1.8, shape = 15, position = position_dodge(width = .6)) +
  geom_errorbarh(height= 0, position = position_dodge(width = .6)) +
  theme_bw(base_size = 8) +
  scale_colour_manual(values = scale_colors_fire) +
  xlim(-6, 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), legend.position = "none", 
        axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold")) +
  labs(title = "Heat shock treatments", x = "Effect size") +
  facet_wrap(~Variable)

figure6b

# Now, we will prepare the figure for the heat data.

# First, we will filter the database.

figure_meta_fire_heat <- figure_meta_fire %>%
  filter(Cue == "Heat")

# Ensure that the moderatos appear in the order we want.

figure_meta_fire_heat$Mods <- factor(figure_meta_fire_heat$Mods, 
                                          levels = c("Seed mass",
                                                     "Non-dormant",
                                                     "Dormant",
                                                     "Dormancy",
                                                     "Xeric",
                                                     "Mesic/Xeric",
                                                     "Microhabitat",
                                                     "Widespread",
                                                     "Restricted",
                                                     "Distribution",
                                                     "Shrub",
                                                     "Herb",
                                                     "Growth form",
                                                     "All"))

figure_meta_fire_heat$Variable <- factor(figure_meta_fire_heat$Variable, 
                                         levels = c("Percentage",
                                                    "Time"))

figure6c <- ggplot(data = figure_meta_fire_heat, aes(y = Mods,
                                                     x = post_mean,
                                                     xmin = Lower,
                                                     xmax = Upper,
                                                     col = Significant)) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.8) +
  geom_point(size = 1.8, shape = 15, position = position_dodge(width = .6)) +
  geom_errorbarh(height= 0, position = position_dodge(width = .6)) +
  theme_bw(base_size = 8) +
  scale_colour_manual(values = scale_colors_fire) +
  xlim(-6, 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), legend.position = "none", 
        axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold")) +
  labs(title = "Heat", x = "Effect size") +
  facet_wrap(~Variable)

figure6c

# Finally, we will prepare the figure with the effect size for smoke treatments.

# First, we will filter the database.

figure_meta_fire_smoke <- figure_meta_fire %>%
  filter(Cue == "Smoke")

figure_meta_fire_smoke$Mods <- factor(figure_meta_fire_smoke$Mods, 
                                           levels = c("Seed mass",
                                                      "Non-dormant",
                                                      "Dormant",
                                                      "Dormancy",
                                                      "Xeric",
                                                      "Mesic/Xeric",
                                                      "Microhabitat",
                                                      "Widespread",
                                                      "Restricted",
                                                      "Distribution",
                                                      "Shrub",
                                                      "Herb",
                                                      "Growth form",
                                                      "All"))

figure_meta_fire_smoke$Variable <- factor(figure_meta_fire_smoke$Variable, 
                                          levels = c("Percentage",
                                                     "Time"))

figure6d <- ggplot(data = figure_meta_fire_smoke, aes(y = Mods,
                                                      x = post_mean,
                                                      xmin = Lower,
                                                      xmax = Upper,
                                                      col = Significant)) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.8) +
  geom_point(size = 1.8, shape = 15, position = position_dodge(width = .6)) +
  geom_errorbarh(height= 0, position = position_dodge(width = .6)) +
  theme_bw(base_size = 8) +
  scale_colour_manual(values = scale_colors_fire) +
  xlim(-6, 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), legend.position = "none", 
        axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold")) +
  labs(title = "Smoke", x = "Effect size") +
  facet_wrap(~Variable)

figure6d

figure6 <- plot_grid(figure6a, figure6b, figure6c, figure6d, nrow = 2, ncol = 2,
                     labels = "AUTO")

figure6

ggsave(filename = "Figure6.svg", device = "svg", height = 12, width = 15, 
       units = "cm", dpi = 600)
