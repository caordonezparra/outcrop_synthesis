# This code was designed to create the figures of our paper: Ordóñez-Parra et al. Seed functional ecology in Brazilian rock outcrop vegetation: an integrative synthesis.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: February 8th, 2024.

# Code run with R v. 4.3.2 on RStudio v. 2023.09.1

#### 0. Set working directory ####

setwd("/Volumes/Apollo M100/2. Doutorado/Projects/outcrop_synthesis")

#### 0. Loading packages (assuming they are already installed) ####

library(readr) # v. 2.1.5 for loading datasets.
library(dplyr) # v. 1.1.4 for data frame manipulation.
library(tidyr) # v. 1.3.0 for cleaning up the dataset.
library(venn) # v. 1.12. for creating Venn diagrams.
library(ggplot2) # v. 3.4.4 for creating different visualizations
library(svglite) # v. 2.1.3 for exporting figures as .svg
library(cowplot) # v. 1.1.2 for arranging plot in grids
library(ggnewscale) # v. 0.4.9 for adding additional axis in graphs.
library(BiocManager) # v. 3.18 for accessing different Bioinformatic packages.
library(treedataverse) # v. 0.0.1 a Bioconductor metapackage with different functions to visualize phylogenetic trees.
library(lcvplants)
library(LCVP)

#### 1. Loading datasets ####

# To create the figures, we need three datasets. First, we'll need the "vegetation.csv" dataset, which contains information on the number of studies (nStudies) and species (nSpecies) associated to each of the four vegetation types present in the "Rock n' Seeds database (https://doi.org/10.1002/ecy.3852).

vegetation <- read_csv("vegetation.csv")

# First, we are uploading the traits dataset, which contains all of the information
# regarding seed trait information.

traits <- read.csv("traits.csv")

# Second, we are loading the topics dataset. This dataset classifies each study, annotating which aspects was assessed for each species on each study.

topics <- read.csv("topics.csv")

# Finally, we are loading the phylogeny we created for our species.

full_tree <- read.tree(file = "outcrop_phylo.tre")

#### 2. Tidying up the trait dataset ####

#### 3. Filtering the trait database ####

# First, obtaning the species accepted name following LCVP (Freiberg et al. 2020), given that these were the names used in the phylogeny and we need the species names in our trait data to match those in the phylogeny. We are saving this  object as "lcvp_names". This procedure will also allow us to level up subspecies and varieties names to the species level and will aid fitering out taxa identified up to genus.

# We are saving this dataset manually to add the name of Pleroma marumbiensis, manually. This species is not present in the LCVP database, but in Flora and Funga do Brasil (http://floradobrasil.jbrj.gov.br/; Last checked, June 16th,  2022). While the accepted name according to Flora do Brasil would be P. marumbiensis, All Tibouchina and Pleroma species where renamed as Pleroma. Therefore, we kept the name of this species as Tibouchina marumbiensis.

lcvp_names <- lcvp_search(traits$Species_reported) %>%
  dplyr::select(Search, Output.Taxon) %>%
  separate(Output.Taxon, into = c("Species", "Genus"), sep = " ") %>%
  unite(col = Species_acceptedLCVP, Species, Genus, sep = "_") %>%
  rename(Species_reported = Search) %>%
  distinct() %>%
  mutate(Species_acceptedLCVP = case_when(Species_reported == "Pleroma marumbiense" ~ "Tibouchina_marumbiensis",
                                          TRUE ~ Species_acceptedLCVP))

# Now, we are going to filter the trait database just to keep the species name as reported in the papers and the seven traits we are interested in (i.e., dispersal period, dispersal syndrome, percentage of embryoless seeds, percentage of viable seeds, seed dry mass, seed water content and primar dormancy). We are also keeping information on species Growth form, geographic distribution and microhabitat.

traits_filtered <- traits %>%
  dplyr::select(Order, Family, Genus, Species_reported, Growth_form, Distribution, Microhabitat, Dispersal_period, 
                Dispersal_syndrome, Embryoless, Viable, Dry_mass, Water_content, Dormancy)

# Now, we join this filtered trait dataset with the table with the accepted names following the LCVP.

traits_joined <- traits_filtered %>%
  left_join(lcvp_names_checked)

#### 2. Figure 1e and 1f ####

# Figure 1 aims to show the different vegetation types included in the study as well as the number of publications and species studied on each of them. Figures A-D correspond to photographs that were arranged manually in Adobe Illustrator. The following code was used to elaborate Figure 1E and F, which are barplots showing the proportion of studies and species from each vegetation type (as a percentage of the total of species and studies).

# We will specify that the column "vegetation" is a factor and create a 'dummy' column that will help us to create a single bar plot.

vegetation <- vegetation %>%
  mutate(Vegetation = as.factor(Vegetation),
         dummy = "dummy")

# Now, we'll indicate the order that we want each vegetation type to appear.

vegetation$Vegetation <- ordered(vegetation$Vegetation,
                                 levels = c("Inselberg",
                                            "Campo de altitude",
                                            "Canga",
                                            "Campo rupestre"))

# We'll also define a vector with the colors we want to use, which were taken directly from the photographs used to illustrate the different vegetation types.

vegetationColors <- c("Inselberg" = "#DECC63", "Campo de altitude" = "#799664",
                      "Canga" = "#997054", "Campo rupestre" = "#B3B4AE")

# Now the dataset is ready, we'll proceed to create the graph with the number of studies

Figure_1E1 <- ggplot(data = figure1e, aes(x = dummy, y = nStudies, fill = Vegetation)) +
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() + 
  theme_classic(base_size = 10) + 
  scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x*100)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(legend.position = "none", axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.title.x = element_text(face = "bold")) +
  ylab("Studies Percentage (%)") +
  scale_fill_manual(values = vegetationColors)

Figure_1E1

# Now, we'll use a similar code to build the same barplot, this time with the percentage of taxa.

Figure_1E2 <- ggplot(data = figure1e, aes(x = dummy, y = nSpecies, fill = Vegetation)) +
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() + 
  theme_classic(base_size = 10) + 
  scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x*100)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(legend.position = "bottom", axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.title.x = element_text(face = "bold")) +
  ylab("Taxa Percentage (%)") +
  scale_fill_manual(values = vegetationColors)

Figure_1E2

# Now, we will arrange these two plots in a grid.

Figure_1E <- plot_grid(Figure_1E1, Figure_1E2, # The plots we want to put in the grid.
                       ncol = 1, nrow = 2, # The number of columns and rows we want in our grid, respectively.
                       labels = c("E", "F"), # The labels of our plots. 
                       rel_heights = c(1,1.75)) # The relative height of each plot. Since Figure 1F has the legend, we need to adjust its height so that bars have the same height.

Figure_1EF

# Now that we're satisfied with the plot, we're going to save this result as a .svg to add it to the photogrids.

ggsave("Figure1EF.svg", width = 17.25, height = 5.5, units = "cm", dpi = 600)

#### 3. Figure 2 ####

# Figure 2 is a phylogenetic tree annotated with the seven seed traits whose phylogenetic signal and structure was assessed. The inner circles represent the different states for our categorical traits: presence/absence of primary dormancy, dispersal syndrome and dispersal season. The outer circles with the bars, illustrate the variation in the quantitative traits: seed dry mass, seed water content, percentage of embryoless seeds and percentage of viable seeds.

# The first step for this figure is to create a database for each trait.

df1 <- traits_joined %>% # Dormancy
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
                                       sum == 1 & LR == 1 ~ "LR"))

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

Figure2 <- ggtree(tr = full_tree, layout = "fan", open.angle = 10) + 
  new_scale_fill() + geom_fruit(data = df1, geom = geom_tile,
                                mapping = aes(y=Species_acceptedLCVP, fill = Dormancy_PA),
                                pwidth= 10, offset = 0.05) + 
  scale_fill_manual(values = c("#EF7E87", "#A95DA2","#D26691")) +
  new_scale_fill() + geom_fruit(data = df2, geom = geom_tile, 
                                mapping = aes(y=Species_acceptedLCVP, fill = Dispersal_syndrome),
                                pwidth= 10, offset = 0.06) + 
  scale_fill_manual(values = c("#FDC300", "#6FAE58","#EA4743")) +
  new_scale_fill() +  geom_fruit(data = df3, geom = geom_tile,
                                 mapping = aes(y = Species_acceptedLCVP, fill = Dispersal_period2),
                                 pwidth = 8, offset = 0.06) +
  scale_fill_manual(values = c("#FB8300", "#28A6A5", "#FFA849", "#53D5D4", "#F2F2F2"))  +
  new_scale_fill() + geom_fruit(data = df4, geom = geom_bar,
                                mapping = aes(y = Species_acceptedLCVP, x = log(Mass + 1)),
                                stat = "identity", offset = 0.06,
                                axis.params = list(axis       = "x",
                                                   text.size  = 1.5,
                                                   hjust      = 0.5,
                                                   vjust      = 1.5,
                                                   nbreak     = 3),
                                grid.params = list(), fill = "#A6A6A6") +
  new_scale_fill() + geom_fruit(data = df5, geom = geom_bar,
                                mapping = aes(y = Species_acceptedLCVP, x = Viable),
                                stat = "identity", offset = 0.06,
                                axis.params = list(axis       = "x",
                                                   text.size  = 1.5,
                                                   hjust      = 0.5,
                                                   vjust      = 1.5,
                                                   nbreak     = 3),
                                grid.params = list(), fill = "#A6A6A6") +
  new_scale_fill() + geom_fruit(data = df6, geom = geom_bar,
                                mapping = aes(y = Species_acceptedLCVP, x = Embryoless),
                                stat = "identity", offset = 0.06,
                                axis.params = list(axis       = "x",
                                                   text.size  = 1.5,
                                                   hjust      = 0.5,
                                                   vjust      = 1.5,
                                                   nbreak     = 3),
                                grid.params = list(), fill = "#A6A6A6") +
  new_scale_fill() + geom_fruit(data = df7, geom = geom_bar,
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

ggsave(filename = "Figure2.svg", device = "svg", dpi = 600)
