# The following codes was designed to create the object "traits_joined", a dataset that contains the trait information data as well as the standardized taxonomy following the Leipzig Catalogue of Vascular Plants (LCVP) (Freiberg et al. 2020) (https://doi.org/10.1038/s41597-020-00702-z).

# Author: Carlos A. Ordóñez Parra (carlos.ordonez.parra@gmail.com)
# Last update: March 4th, 2024.

# This code was run using R v. 4.3.2 in R Studio v. 2023.09.1

#### 0. Set working directory ####

setwd("/Volumes/Apollo M100/2. Doutorado/Projects/outcrop_synthesis/Preprocessing")

#### 1. Load the packages (assuming they are already installed) ####

library(readr) # v. 2.1.5 for loading datasets.
library(dplyr) # v. 1.1.4 for data frame manipulation.
library(tidyr) # v. 1.3.0 for cleaning up the dataset.
library(lcvplants) # 2.1.0 for assessing the LCVP
library(LCVP) # 3.0.1 for functions to check species name against the LCVP

#### 2. Load the dataset ####

# For this code, we only need the trait dataset, which corresponds to the trait.csv file in Ordóñez-Parra et al. (2023) (https://doi.org/10.1002/ecy.3852)

traits <- read_csv("/Volumes/Apollo M100/2. Doutorado/Projects/outcrop_synthesis/traits.csv")

#### 3. Tidying up the trait dataset ####

# First, we'll obtaning the species accepted name following LCVP, given that these were the names used in the phylogeny (see 'phylogeny.R') and we need the species names in our trait data to match those in the phylogeny. We are saving this  object as "lcvp_names". This procedure will also allow us to level up subspecies and varieties names to the species level and will aid fitering out taxa identified up to genus.

# We are saving this dataset manually to add the name of Pleroma marumbiensis, manually. This species is not present in the LCVP database, but it is reported in Flora and Funga do Brasil (http://floradobrasil.jbrj.gov.br/; Last checked, June 16th,  2022). While the accepted name according to Flora do Brasil would be Pleroma marumbiensis, All Tibouchina and Pleroma species where renamed as Pleroma. Therefore, we kept the name of this species as Tibouchina marumbiensis.

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
  left_join(lcvp_names) %>%
  filter(Species_acceptedLCVP != "NA_NA")

# We are saving this file as a .csv in case we need to load it in future analysis.

write.csv(traits_joined, "traits_joined.csv")
