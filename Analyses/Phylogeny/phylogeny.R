# This code was designed to update the species nomenclature following the Leipzig
# Catalogue for Vascular Plants (Freiberg et al. 2020) and create the phylogeny 
# required for the analysis of this paper.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)

# Last updated in November 7th, 2022.

# Code run with R version 4.2.0. in RStudio

#### 0. Loading packages #####

# readr v. 2.1.3. for loading the dataset

if(!require(readr)){
  install.packages("readr")
  require(readr)
}

# dplyr v. 1.0.10 for dataset manipulation

if(!require(dplyr)){
  install.packages("dplyr")
  require(dplyr)
}

# tidyr v. 1.2.1 for cleaning up the dataset.

if(!require(tidyr)){
  install.packages("tidyr")
  require(tidyr)
}

# The packages LCVP, lcvpplants and V.Phylomaker2 were installed via GitHub using
# the following code (Last run, June 16th 2022).

if(!require(devtools)){
  install.packages("devtools")
  require(devtools)
}

install_github("idiv-biodiversity/LCVP")
install_github("idiv-biodiversity/lcvplants")
install_github("jinyizju/V.PhyloMaker2")

# LCVP v. 2.0. to access the Leipzig Catalogue of Vascular Plants (LCVP)

library(LCVP)

# lcvplants v. 2.0.0 for the tools to update nomenclature following LCVP

library(lcvplants)

# V.PhyloMaker2 v. 0.1.0 for gerating the phylogeny.

library(V.PhyloMaker2)

#### 1. Loading datasets ####

# We are obtaning our species name frome the "traits.csv" from the Rock n' Seeds
# database (Ordóñez-Parra et al. 2022).

traits <- read.csv("traits.csv")

#### 2. Generating the species list ####

# First, we are using the lcvp_search function from the lcvplants to check for the
# updated names according to LCVP

specieslcvp <- lcvp_search(traits$Species_reported)

# Since several species appear more than once in our dataset, we are keeping only
# distinct observations. Also, we filtering out taxa identified only to genera level. 
# Also, the name Pleroma marumbiense is not present in LCVP. Yet it appears as an 
# accepted name in Flora do Brasil (http://floradobrasil.jbrj.gov.br/; Last checked,
# June 16th, 2022). Therefore, we are inserting it manually on the "Output.Taxon" column.

specieslcvp2 <- specieslcvp %>%
  dplyr::select(Output.Taxon, Family) %>%
  distinct() %>%
  add_row(Output.Taxon = "Tibouchina marumbiense", Family = "Melastomataceae")

# Now, in order to create the file we need to create the phylogeny using V.PhyloMaker2
# we are separating the species name into Genus and Species, and getting rid of
# the authors names.

specieslcvp3 <- specieslcvp2 %>%
  separate(col = Output.Taxon, into = c("Genus", "Species"), sep = " ", 
           extra = "drop" ) %>%
  unite(col = "Genus_Species", Genus, Species,sep = "_", remove = FALSE) %>%
  dplyr::select(-Species) %>%
  rename(Species = Genus_Species) %>%
  distinct()

# Now, we have the final database to use for building the phylogeny. Still, preliminar
# attempts of creating the phylogeny showed that several of the species in our
# list did not have a close relative in the phylogeny neither their genera was
# present. Therefore, we are filtering out this species to later add them as 
# polytomies in their corresponding position.

filtered_species <- c("Austrocritonia_velutina", "Cavalcantia_glomerata", 
                      "Cavalcantia_percymosa", "Parapiqueria_cavalcantei")

specieslcvp4 <- specieslcvp3 %>%
  filter(!Species %in% filtered_species,
         Species  != "NA_NA") %>%
  distinct()

# We are saving this final version of the species list to add the columns 
# species.relative and genus.relative to bind two of our species with a closely
# related species that is present in the phylogeny. Specifically, we are  binding
# Brasilianthus carajensis and Monogereiron carajensis with their relatives Nepsera 
# aquatica and Heterocondylus pumilus, respectively (see Almeda et al. 2016, Rivera
# et al. 2016). We are also chaning the family of Myrsine gardeniana to Primulaceae,
# since that is the name present in the megaphylogeny.

specieslcvp_final <- specieslcvp4 %>%
  mutate(species.relative = case_when(Species == "Brasilianthus_carajensis" ~ "Nepsera_aquatica",
                                      Species == "Monogereion_carajensis" ~ "Heterocondylus_pumilus",
                                      TRUE ~ "NA"),
         genus.relative = case_when(Species == "Brasilianthus_carajensis" ~ "Nepsera",
                                    Species == "Monogereion_carajensis" ~ "Heterocondylus",
                                    TRUE ~ "NA"),
         Family = case_when(Family == "Myrsinaceae" ~ "Primulaceae",
                            TRUE ~ Family)) %>%
  rename("species" = "Species",
         "genus" = "Genus",
         "family" = "Family")

#### 3. Generating the phylogeny ####

# Now we proceed to generate the phylogeny. First, we bind B. carajensis and M.
# carajensis with their respective relatives.

binded <- bind.relative(specieslcvp_final, tree = GBOTB.extended.LCVP, 
                        nodes = nodes.info.1.LCVP)

# Then, we use this phylogeny that has this species binded and use it to generate
# our phylogeny by keeping only the species in our list.

phylo <- phylo.maker(sp.list = binded$species.list, tree = binded$phylo,
                     nodes = binded$nodes.info, scenarios = "S3")

plot.phylo(phylo$scenario.3, type = "fan", cex = 0.3, label.offset = 1.5, 
           no.margin = TRUE)

# Now, we are proceeding to insert those species we filtered at the begining as
# polytomies. They are all from the Eupatorinae tribe in Asteraceae. In this tribe
# relationships between genera are unresolved. Therefore, we are including them
# as polytomies in the node of the clade formed by Mikania, Praxelis, Chromolaena,
# Monogereiron and Heterocondylus.

write.tree(phylo$scenario.3, "phylo1.tre")

phylo2 <- at.node(phylo$scenario.3, location.node = phylo$scenario.3$edge[17], tip.label = "Austrocritonia_velutina")
phylo3 <- at.node(phylo2, location.node = phylo$scenario.3$edge[18], tip.label = "Cavalcantia_glomerata")
phylo4 <- at.node(phylo3, location.node = phylo$scenario.3$edge[19], tip.label = "Parapiqueria_cavalcantei")
phylo5 <- at.node(phylo4, location.node = phylo$scenario.3$edge[22], tip.label = "Cavalcantia_percymosa")

phylo_figure <- plot.phylo(phylo5, type = "fan", cex = 0.3, label.offset = 1.5, 
                           no.margin = TRUE)

write.tree(phylo5, "outcrop_phylo.tre")

# Now, we have the final phylogeny we are going to use for our analyses.

#### 4. References #####

# Almeda et al. (2016) https://doi.org/10.11646/phytotaxa.273.4.3
# Freiberg et al. (2020) https://doi.org/10.1038/s41597-020-00702-z
# Rivera et al. (2016) https://doi.org/10.1016/j.ympev.2015.11.013
