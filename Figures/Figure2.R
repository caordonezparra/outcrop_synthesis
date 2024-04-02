# The following code is designed to create Figure 2 of Ordóñez-Parra et al, which is a is a phylogenetic tree annotated with the seven seed traits whose phylogenetic signal and was assessed. 

# The inner circles represent the different states for our three categorical traits: presence/absence of primary dormancy, dispersal syndrome and dispersal season. The outer circles with the bars, illustrate the variation in the quantitative traits: seed dry mass, seed water content, percentage of embryoless seeds and percentage of viable seeds.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: April 2nd, 2024.

#### 0. Load the packages (assuming they are already installed) ####

library(gert) # v. 2.0.1 for managing Git repositories
library(readr) # v. 2.1.5 for loading datasets.
library(ape) # v. 5.7-1 for loading the phylogenetic tree.
library(dplyr) # v. 1.1.4 for data frame manipulation.
library(tidyr) # v. 1.3.0 for cleaning up the dataset.
library(ggplot2) # v. 3.4.4 for creating different visualizations
library(svglite) # v. 2.1.3 for exporting figures as .svg
library(ggnewscale) # v. 0.4.10 for adding additional axis in graphs.
library(BiocManager) # v. 3.18 for accessing different Bioinformatic packages.
library(treedataverse) # v. 0.0.1 a Bioconductor metapackage with different functions to visualize phylogenetic trees

##### 1. Set working directory ####

setwd(str_c(git_find(),"/Figures"))

#### 2. Load the dataset and the phylogeny ####

# For creating this figure, we will be using the traits_joined data set, which contains both the information of the seed functional traits and the updated taxonomy following the Leipzig Catalogue of Vascular Plants. (See Section 2 in phylogenetic_signal.R)

traits_joined <- read_csv(str_c(git_find(),"/Preprocessing/traits_joined.csv"))

# We will also need the phylogenetic tree we have created in phylogeny.R, which is available in "outcrop_phylo.tre"

full_tree <- read.tree(file = str_c(git_find(), "/outcrop_phylo.tre"))

#### 3. Prepare the datasets for the Figure ####

# The first step for creating this figure is to create a datasets with the mean value of each trait for each species. With seed dormancy and dispersal period, we'll be doing some pre-processing to plot the cases where a species presents record for D and ND or for more than one dispersal period.

# Primary dormancy

df1 <- traits_joined %>% # Take the 'trait_joined' data set.
  select(Species_acceptedLCVP, Dormancy) %>% # Select these two columns.
  filter(!is.na(Species_acceptedLCVP), # Filter out species with no species record.
         !is.na(Dormancy), # Filter out species without information on primary dormancy.
         Dormancy != "NC") %>% # Filter out non-conclusive 'NC' dormancy records.
  distinct() %>% # Keep only distinct observations.
  mutate(Record = 1) %>% # Create a new column called 'Record' where every observation has 1 as value.
  pivot_wider(names_from = Dormancy, values_from = Record, values_fill = 0) %>% # Reshape the data set to take the     Dormancy column and create two new columns (based on 'Dormancy'). Fill each observation with the value of 'Rercord'   or zero if no 'Record' value is available.
  mutate(Dormancy_PA = case_when(ND == 1 & D == 0 ~ "ND", # Create a new column called 'Dormancy_PA' based on the      values of the columns we have just created.
                                 ND == 0 & D == 1 ~ "D",
                                 ND == 1 & D == 1 ~ "Both"))

# Dispersal syndrome

df2 <- syndrome <- traits_joined %>% # Take the 'trait_joined' data set.
  select(Species_acceptedLCVP, Dispersal_syndrome) %>% # Select these two columns.
  filter(!is.na(Species_acceptedLCVP), # Filter out species with no species record.
         !is.na(Dispersal_syndrome)) %>% # Filter out species without information on dispersal syndrome.
  distinct() # Keep only distinct observations.

# Dispersal period

df3 <- traits_joined %>% # Take the 'trait_joined' data set.
  select(Species_acceptedLCVP, Dispersal_period) %>% # Select these two columns.
  filter(!is.na(Species_acceptedLCVP), # Filter out species with no species record.
         !is.na(Dispersal_period)) %>% # Filter out species without information on dispersal period.
  distinct() %>% # Keep only distinct observations.
  mutate(Record = 1) %>% # Create a new column called 'Record' where every observation has 1 as value.
  pivot_wider(names_from = Dispersal_period, values_from = Record, values_fill = 0) %>% # Reshape the data set to     take the 'Dispersal_period' column and create four new columns (based on 'Dispersal_period'). It will fill each      observation with the value of 'Rercord' or zero if no 'Record' value is available.
  mutate(sum = ED + ER + LD + LR, # Create a new column with the sum of the values of the four columns we have just   created.
         Dispersal_period2 = case_when(sum > 1 ~ "More than 1", # Create a new column, called 'Dispersal_period'              based on the values on the 'sum' column. This will allow us to identify whether a species disperse their             seeds in more than one period.
                                       sum == 1 & LD == 1 ~ "LD",
                                       sum == 1 & ED == 1 ~ "ED",
                                       sum == 1 & ER == 1 ~ "ER",
                                       sum == 1 & LR == 1 ~ "LR"),
         Dispersal_period2 = as.factor(Dispersal_period2)) # Ensure that R understands this new column as a factor.

# Seed dry mass
  
df4 <- traits_joined %>% # Take the 'trait_joined' data set.
  select(Species_acceptedLCVP, Dry_mass) %>% # Select these two columns.
  filter(!is.na(Species_acceptedLCVP), # Filter out species with no species record.
         !is.na(Dry_mass)) %>% # Filter out species without information on seed dry mass.
  group_by(Species_acceptedLCVP) %>% # Group the subsequent calculations by species.
  summarise(Mass = mean(Dry_mass)) # Calculate the mean seed dry mass.

# Percentage of viable seeds

df5 <- traits_joined %>% # Take the 'trait_joined' data set.
  select(Species_acceptedLCVP, Viable) %>% # Select these two columns.
  filter(!is.na(Species_acceptedLCVP), # Filter out species with no species record.
         !is.na(Viable)) %>% # Filter out species without information on the percentage of viable seeds.
  group_by(Species_acceptedLCVP) %>% # Group the subsequent calculations by species.
  summarise(Viable = mean(Viable)) # Calculate the mean percentage of viable seeds.

# Percentage of embryoless seeds

df6 <- traits_joined %>% # Take the 'trait_joined' data set.
  select(Species_acceptedLCVP, Embryoless) %>% # Select these two columns.
  filter(!is.na(Species_acceptedLCVP), # Filter out species with no species record.
         !is.na(Embryoless)) %>% # Filter out species without information on the percentage of embryoless seeds.
  group_by(Species_acceptedLCVP) %>% # Group the subsequent calculations by species.
  summarise(Embryoless = mean(Embryoless)) # Calculate the mean percentage of viable seeds.

# Seed water content

df7 <- traits_joined %>% # Take the 'trait_joined' data set.
  select(Species_acceptedLCVP, Water_content) %>% # Select these two columns.
  filter(!is.na(Species_acceptedLCVP), # Filter out species with no species record.
         !is.na(Water_content)) %>% # Filter out species without information on seed water content.
  group_by(Species_acceptedLCVP) %>% # Group the subsequent calculations by species.
  summarise(Water_content = mean(Water_content)) # Calculate the mean seed water content.

##### 4. Create the figure ####

# Now that we have created all the datasets, we can proceed to the creating the figure. Most of the following lines is heavily based on the following link: https://yulab-smu.top/treedata-book/chapter7.html

Figure2 <- ggtree(tr = full_tree, layout = "fan", open.angle = 10) + # We specify the phylogenetic tree that is   going to be used as a base (tr argument), how we want to plot the phylogeny (layout = 'fan', which is a circular     and rooted), and the opening angle that will allow us to enter the legends of the chart's axes (open.angle).
  geom_fruit(data = df1, geom = geom_tile, # For categorical traits, we'll use geom_tile to create squares...
             mapping = aes(y = Species_acceptedLCVP, fill = Dormancy_PA), # that will change of color according to     the different levels inside the trait.
             pwidth = 10, offset = 0.05) + # We'll define the size of the square (pwidth) and the separation between   this and the previous geometry.
  scale_fill_manual("Dormancy_PA", values = c("#EF7E87", "#A95DA2", "#D26691")) + # We'll also define manually the     color palette.
  new_scale_fill() + # We'll use this command everytime we want to add a new layer to our graph.
  geom_fruit(data = df2, geom = geom_tile, 
             mapping = aes(y = Species_acceptedLCVP, fill = Dispersal_syndrome),
             pwidth = 10, offset = 0.06) + 
  scale_fill_manual("Dispersal_syndrome", values = c("#FDC300", "#6FAE58", "#EA4743")) +
  new_scale_fill() +
  
  geom_fruit(data = df3, geom = geom_tile,
             mapping = aes(y = Species_acceptedLCVP, fill = Dispersal_period2),
             pwidth = 10, offset = 0.07) +
  scale_fill_manual("Dispersal_period2", values = c("#FB8300", "#28A6A5", "#FFA849", "#53D5D4", "#F2F2F2"))  +
  geom_fruit(data = df4, geom = geom_bar, # For quantitative traits we'll use bar charts, instead.
             mapping = aes(y = Species_acceptedLCVP, x = log(Mass + 1)),
             stat = "identity", offset = 0.06,
             axis.params = list(axis       = "x",
                                text.size  = 1.5,
                                hjust      = 0.5,
                                vjust      = 1.5,
                                nbreak     = 3,
                                limits = c(0, 6)),
             grid.params = list(), fill = "#A6A6A6") +
  geom_fruit(data = df5, geom = geom_bar,
             mapping = aes(y = Species_acceptedLCVP, x = Viable),
             stat = "identity", offset = 0.06,
             axis.params = list(axis       = "x",
                                text.size  = 1.5,
                                hjust      = 0.5,
                                vjust      = 1.5,
                                nbreak     = 2,
                                limits = c(0,100)),
             grid.params = list(), fill = "#A6A6A6") +
  geom_fruit(data = df6, geom = geom_bar,
             mapping = aes(y = Species_acceptedLCVP, x = Embryoless),
             stat = "identity", offset = 0.06,
             axis.params = list(axis       = "x",
                                text.size  = 1.5,
                                hjust      = 0.5,
                                vjust      = 1.5,
                                nbreak     = 2,
                                limits = c(0,100)),
             grid.params = list(), fill = "#A6A6A6") +
  geom_fruit(data = df7, geom = geom_bar,
             mapping = aes(y = Species_acceptedLCVP, x = Water_content),
             stat = "identity", offset = 0.06,
             axis.params = list(axis       = "x",
                                text.size  = 1.5,
                                hjust      = 0.5,
                                vjust      = 1.5,
                                nbreak     = 3,
                                limits = c(0, 30)),
             grid.params = list(), fill = "#A6A6A6") +
  theme(legend.position = "none") # We will eliminate the legend for this graph and we'll add it manually on         Illustrator.

Figure2

# Once we are satisfied with the final result, we can save it as a .svg file to add some details in Illustrator.

ggsave(filename = "Figure2.svg", device = "svg", dpi = 600,
       height = 16.5, width = 16.5, units = "cm")
