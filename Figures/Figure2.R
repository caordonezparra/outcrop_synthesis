# The following code is designed to create Figure 2 of Ordóñez-Parra et al, which is a is a phylogenetic tree annotated with the seven seed traits whose phylogenetic signal and structure was assessed. 

#The inner circles represent the different states for our categorical traits: presence/absence of primary dormancy, dispersal syndrome and dispersal season. The outer circles with the bars, illustrate the variation in the quantitative traits: seed dry mass, seed water content, percentage of embryoless seeds and percentage of viable seeds.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: March 5, 2024.

##### 0. Set working directory ####

setwd("/Volumes/Apollo M100/2. Doutorado/Projects/outcrop_synthesis/Figures")


#### 1. Load the packages (assuming they are already installed) ####

library(readr) # v. 2.1.5 for loading datasets.
library(ape) # v. 5.7-1 for loading the phylogenetic tree.
library(dplyr) # v. 1.1.4 for data frame manipulation.
library(tidyr) # v. 1.3.0 for cleaning up the dataset.
library(ggplot2) # v. 3.4.4 for creating different visualizations
library(svglite) # v. 2.1.3 for exporting figures as .svg
library(ggnewscale) # v. 0.4.10 for adding additional axis in graphs.
library(BiocManager) # v. 3.18 for accessing different Bioinformatic packages.
library(treedataverse) # v. 0.0.1 a Bioconductor metapackage with different functions to visualize phylogenetic trees.

#### 2. Load the dataset and the phylogeny ####

# For creating this figure, we will be using the traits_joined dataset, which contains both the information of the seed functional traits and the updated taxonomy following the Leipzig Catalogue of Vascular Plants. See the traits_joined.R in the Preprocessing folder to see the procedure for defining and producing this dataset.

read_csv("/Volumes/Apollo M100/2. Doutorado/Projects/outcrop_synthesis/Preprocessing/traits_joined.csv")

# We will also need the phylogenetic tree we have created in phylogeny.R, which is available in "outcrop_phylo.tre"

full_tree <- read.tree(file = "/Volumes/Apollo M100/2. Doutorado/Projects/outcrop_synthesis/outcrop_phylo.tre")

#### 3. Making the figure ####

# The first step for creating this figure is to create a database for each trait. In some cases, with seed dormancy and dispersal period, we'll be doing some pre-processing to plot the cases where a species presents record for D and ND or for more than one dispersal period.

df1 <- traits_joined %>% # Primary dormancy
  select(Species_acceptedLCVP, Dormancy) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Dormancy),
         Dormancy != "NC") %>%
  distinct() %>% 
  mutate(Record = 1) %>% 
  pivot_wider(names_from = Dormancy, values_from = Record, values_fill = 0) %>%
  mutate(Dormancy_PA = case_when(ND == 1 & D == 0 ~ "ND",
                                 ND == 0 & D == 1 ~ "D",
                                 ND == 1 & D == 1 ~ "Both"))

df2 <- syndrome <- traits_joined %>% # Dispersal syndrome
  select(Species_acceptedLCVP, Dispersal_syndrome) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Dispersal_syndrome)) %>% 
  distinct()

df3 <- traits_joined %>% # Dispersal_period
  select(Species_acceptedLCVP, Dispersal_period) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Dispersal_period)) %>%
  distinct() %>%
  mutate(Record = 1) %>%
  pivot_wider(names_from = Dispersal_period, values_from = Record, values_fill = 0) %>%
  mutate(sum = ED + ER + LD + LR,
         Dispersal_period2 = case_when(sum > 1 ~ "More than 1",
                                       sum == 1 & LD == 1 ~ "LD",
                                       sum == 1 & ED == 1 ~ "ED",
                                       sum == 1 & ER == 1 ~ "ER",
                                       sum == 1 & LR == 1 ~ "LR"),
         Dispersal_period2 = as.factor(Dispersal_period2))

df4 <- traits_joined %>% # Dry mass
  select(Species_acceptedLCVP, Dry_mass) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Dry_mass)) %>%
  group_by(Species_acceptedLCVP) %>%
  summarise(Mass = mean(Dry_mass))

df5 <- traits_joined %>% # Viable
  select(Species_acceptedLCVP, Viable) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Viable)) %>%
  group_by(Species_acceptedLCVP) %>%
  summarise(Viable = mean(Viable/100))

df6 <- traits_joined %>% #Embryoless
  select(Species_acceptedLCVP, Embryoless) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Embryoless)) %>%
  group_by(Species_acceptedLCVP) %>%
  summarise(Embryoless = mean(Embryoless/100))

df7 <- traits_joined %>% # Water_content
  select(Species_acceptedLCVP, Water_content) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Water_content)) %>%
  group_by(Species_acceptedLCVP) %>%
  summarise(Water_content = mean(Water_content/100))

# Once all these datasets are created, we can proceed with creating the figure.

Figure2 <- ggtree(tr = full_tree, layout = "fan", open.angle = 10) + 
  geom_fruit(data = df1, geom = geom_tile,
             mapping = aes(y = Species_acceptedLCVP, fill = Dormancy_PA),
             pwidth = 10, offset = 0.05) + 
  scale_fill_manual("Dormancy_PA", values = c("#EF7E87", "#A95DA2", "#D26691")) +
  new_scale_fill() +
  
  geom_fruit(data = df2, geom = geom_tile, 
             mapping = aes(y = Species_acceptedLCVP, fill = Dispersal_syndrome),
             pwidth = 10, offset = 0.06) + 
  scale_fill_manual("Dispersal_syndrome", values = c("#FDC300", "#6FAE58", "#EA4743")) +
  # new_scale_fill() +
  
  # geom_fruit(data = df3, geom = geom_tile,
             #mapping = aes(y = Species_acceptedLCVP, fill = Dispersal_period2),
             #pwidth = 10, offset = 0.06) +
  # scale_fill_manual("Dispersal_period2", values = c("#FB8300", "#28A6A5", "#FFA849", "#53D5D4", "#F2F2F2"))  +
  
  geom_fruit(data = df4, geom = geom_bar,
             mapping = aes(y = Species_acceptedLCVP, x = log(Mass + 1)),
             stat = "identity", offset = 0.06,
             axis.params = list(axis       = "x",
                                text.size  = 1.5,
                                hjust      = 0.5,
                                vjust      = 1.5,
                                nbreak     = 3),
             grid.params = list(), fill = "#A6A6A6") +
  geom_fruit(data = df5, geom = geom_bar,
             mapping = aes(y = Species_acceptedLCVP, x = Viable),
             stat = "identity", offset = 0.06,
             axis.params = list(axis       = "x",
                                text.size  = 1.5,
                                hjust      = 0.5,
                                vjust      = 1.5,
                                nbreak     = 3),
             grid.params = list(), fill = "#A6A6A6") +
  geom_fruit(data = df6, geom = geom_bar,
             mapping = aes(y = Species_acceptedLCVP, x = Embryoless),
             stat = "identity", offset = 0.06,
             axis.params = list(axis       = "x",
                                text.size  = 1.5,
                                hjust      = 0.5,
                                vjust      = 1.5,
                                nbreak     = 3),
             grid.params = list(), fill = "#A6A6A6") +
  geom_fruit(data = df7, geom = geom_bar,
             mapping = aes(y = Species_acceptedLCVP, x = Water_content),
             stat = "identity", offset = 0.06,
             axis.params = list(axis       = "x",
                                text.size  = 1.5,
                                hjust      = 0.5,
                                vjust      = 1.5,
                                nbreak     = 3),
             grid.params = list(), fill = "#A6A6A6") +
  theme(legend.position = "bottom")

Figure2

# Once we are satisfied with the final result, we can save it as a .svg file.

ggsave(filename = "Figure2.svg", device = "svg", dpi = 600)