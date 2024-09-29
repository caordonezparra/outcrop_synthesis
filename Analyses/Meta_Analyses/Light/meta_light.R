# The following code was designed to assess the effect of light on germination percentage by performing meta-analyses through fitting binomial phylogenetic generalized mixed models with Bayesian estimation using Markov chain Monte Carlo (MCMC) method as implemented in the MCMCglmm R-package.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: April 18th, 2024.

#### 0. Loading packages (assuming they are already installed) ####

library(gert) # v. 2.0.1 for managing Git repositories.
library(stringr) # v. 1.5.1 for string manipulation.
library(readr) # v. 2.1.5. for database loading.
library(dplyr) # v. 1.1.4 for database manipulation.
library(tidyr) # v. 1.3.0. (also) for database manipulation
library(tibble) # v. 3.2.1 (also) for database manipulation
library(picante) # v. 1.8.2 for loading phylogenetic trees
library(MCMCglmm) # v. 2.35 for generalized mixed models with Bayesian estimation using MCMC method.

#### 1. Setting working directory ####

setwd(str_c(git_find(),"/Analyses/Meta_Analyses/Light"))

#### 2. Loading datasets ####

# For this script we will be using four data sets. First, we are loading a data frame that classifies the studies into the topics they studied. For this classification three columns are used: Factor, Major and Minor.

topics <- read_csv(str_c(git_find(),"/topics.csv"))

# Second, we are uploading a database with the germination data provided by the authors.

germination <- read_csv(str_c(git_find(),"/germination.csv"))

# Third, we are uploading the database of seed functional traits, which includes information on species distribution, microhabitats, growth form and seed mass.

traits <- read_csv(str_c(git_find(),"/traits.csv"))

# Finally, we are uploading a database with the species names standardized against the Leipzig Catalogue of Vascular Plants (Freiberg et al. 2020)

lcvp_names_checked <- read.csv(str_c(git_find(),"/lcvp_names_checked.csv"))

# We are also uploading the phylogenetic tree with all the species in our dataset.

full_tree <- read.tree(str_c(git_find(),"/Analyses/Phylogeny/outcrop_phylo.tre"))

#### 2. Getting the data from studies that assessed the effect of light ####

# The following code creates a data set of the germination experiments that tested the effect of light availability. The control treatment in all cases were seeds incubated under dark conditions.First we're getting the list of studies researching the effects of light availabilty on germination. Therefore, we are filtering the "topics" database looking for studies where "Major" is "Light" and "Minor" is "Availability".

light <- topics %>%
  filter(Major == "Light",
         Minor == "Availability")

# Then, we are getting the IDs of each of these studies.

light_studies <- unique(light$ID)
light_studies

# Now, we filter the germination database looking for:

light_experiments <- germination %>%
  filter(ID %in% light_studies, # Studies with whose IDs are in "light_studies",
         Chemical_Compound == "Control" | is.na(Chemical_Compound), # Experiments that only tested light effects.
         Storage_Time == "Control" | is.na(Storage_Time),
         Scarification == "Sandpaper" | is.na(Scarification),
         Gut_Passage == "Extracted seeds"| is.na(Gut_Passage),
         HeatShock_Temperature == "Control" | is.na(HeatShock_Temperature),
         Smoke == "Control" | is.na(Smoke)) %>%
  mutate(Condition = case_when(Photoperiod == 0 ~ "Dark",
                               Photoperiod != 0 ~ "Light")) %>%
  dplyr::select(1:3, Condition, Temperature, GermSeeds, nSeeds) %>%
  left_join(traits) %>%
  dplyr::select(1:7, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits we're interested in.
  left_join(lcvp_names_checked) %>%
  filter(Growth_form == "Herb" | Growth_form == "Shrub")

# to avoid the confounding effects from treatments with extremely low germination, we only included observation where either treatment (light and darkness) had >10% germination. For that purpose, we will create a new object with such observations

exclude_light_experiments <- light_experiments %>%
  dplyr::select(1:7) %>%
  mutate(GermProp = GermSeeds/nSeeds) %>%
  group_by(ID, Species_reported, Species_accepted, Condition, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Condition, values_from = meanGermProp) %>%
  filter(is.na(Dark) | Dark < 0.1 & Light < 0.1) %>%
  dplyr::select(1:4)

light_experiments2 <- light_experiments %>%
  anti_join(exclude_light_experiments) %>%
  filter(!(Dormancy == "D" & !(ID %in% c(14, 78))), 
         Dormancy != "NC") # We are keeping studies with ID 14 and 78, since this studies successfuly alliviated           dormancy through scarification.

# Just to have an idea of the size of our dataset, we will look and how many species and studies we have.

n_distinct(light_experiments2$ID) # We have 16 studies
n_distinct(light_experiments2$Species_acceptedLCVP) # For a total of 71 species.

# Also, just to confirm which studies are being include in the meta-analysis, we will create the following object.

included_light_studies <- light_experiments2 %>% # Take the light_experiments2 dataset
  dplyr::select(ID) %>% # Select the 'ID" column
  distinct() # And keep only unique records (i.e., if the ID is repeated, keep only the first line).

included_light_studies

# Before proceeding to the next step, we need to input missing seed mass data. Thus, we'll be exporting the light_experiments2 dataset, to manually add the seed mass values.

write.csv(light_experiments2, file = "light_experiments_missing.csv")

# After imputing the missing seed mass data manually, we'll be loading the final dataset that we are going to use to fit the model. We'll be also adding a new column call n, that will allow us to control the effect of various observation between studies.

light_experiments_final <- read.csv(str_c(git_find(),
                                          "/Analyses/Meta_Analyses/Light/light_experiments_final.csv")) %>%
  rename(Light = Condition,
          animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Light = case_when(Light == "Dark" ~ 0,
                           Light == "Light" ~ 1),
         Light = as.factor(Light))

# We will also need to prune the phylogenetic tree, so that it only keeps the observations present in the 'light_experiments_final' dataset.

meta_light_names <- light_experiments_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_light_names <- names_for_pruning(meta_light_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_light_tree2 <- prune.sample(meta_light_names, phylo = full_tree)
meta_light_tree2$node.label <- NULL

# Now, we will set the number of desired iterations based on De Villemereuil & Nakagawa (2014)

nite = 500000
nthi = 50
nbur = 50000

# Set priors based on https://doi.org/10.1017/S0960258517000332

priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

# First, we are running the MCMCglmm with the whole dataset.

light_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Light, 
                                random = ~ animal + ID + n:ID, 
                                family = "multinomial2", pedigree = meta_light_tree, 
                                prior = priors, data = light_experiments_final,
                                nitt = nite, thin = nthi, burnin = nbur, 
                                verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                pr = FALSE, pl = FALSE)

summary(light_all)
n_distinct(light_experiments_final$animal)

lambda <- light_all$VCV[,'animal']/
  (light_all$VCV[,'animal']+light_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Now, we'll do the analysis using only herbs.

light_experiments_herbs <- light_experiments_final %>%
  filter(Growth_form == "Herb")

light_herbs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Light, 
                                random = ~ animal + ID + n:ID, 
                                family = "multinomial2", pedigree = meta_light_tree, 
                                prior = priors, data = light_experiments_herbs,
                                nitt = nite, thin = nthi, burnin = nbur, 
                                verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                pr = FALSE, pl = FALSE)

summary(light_herbs)
n_distinct(light_experiments_herbs$animal)

lambda <- light_herbs$VCV[,'animal']/
  (light_herbs$VCV[,'animal']+light_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Only shrubs

light_experiments_shrubs <- light_experiments_final %>%
  filter(Growth_form == "Shrub")

light_shurbs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Light, 
                                  random = ~ animal + ID + n:ID, 
                                  family = "multinomial2", pedigree = meta_light_tree, 
                                  prior = priors, data = light_experiments_shrubs,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

summary(light_shurbs)
n_distinct(light_experiments_shrubs$animal)

lambda <- light_shurbs$VCV[,'animal']/
  (light_shurbs$VCV[,'animal']+light_shurbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Only species restricted to the outcrop vegetation.

light_experiments_restricted <- light_experiments_final %>%
  filter(Distribution == "Restricted")

light_restricted <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Light, 
                                   random = ~ animal + ID + n:ID, 
                                   family = "multinomial2", pedigree = meta_light_tree, 
                                   prior = priors, data = light_experiments_restricted,
                                   nitt = nite, thin = nthi, burnin = nbur, 
                                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                   pr = FALSE, pl = FALSE)

summary(light_restricted)
n_distinct(light_experiments_restricted$animal)

lambda <- light_restricted$VCV[,'animal']/
  (light_restricted$VCV[,'animal']+light_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Only species with widespread distribution.

light_experiments_widespread <- light_experiments_final %>%
  filter(Distribution == "Widespread")

light_widespread <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Light, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_light_tree, 
                                       prior = priors, data = light_experiments_widespread,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

summary(light_widespread)
n_distinct(light_experiments_widespread$animal)

lambda <- light_widespread$VCV[,'animal']/
  (light_widespread$VCV[,'animal']+light_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Only species from mesic microhabitats.

light_experiments_mesic <- light_experiments_final %>%
  filter(Microhabitat == "Mesic")

light_mesic <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Light, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_light_tree, 
                                       prior = priors, data = light_experiments_mesic,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

summary(light_mesic)
n_distinct(light_experiments_mesic$animal)

lambda <- light_mesic$VCV[,'animal']/
  (light_mesic$VCV[,'animal']+light_mesic$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Species present only on Xeric environments.

light_experiments_xeric <- light_experiments_final %>%
  filter(Microhabitat == "Xeric")

light_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Light, 
                                  random = ~ animal + ID + n:ID, 
                                  family = "multinomial2", pedigree = meta_light_tree, 
                                  prior = priors, data = light_experiments_xeric,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

summary(light_xeric)
n_distinct(light_experiments_xeric$animal)

lambda <- light_xeric$VCV[,'animal']/
  (light_xeric$VCV[,'animal']+light_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Species present in both kinds of microhabitat.

light_experiments_mesic_xeric <- light_experiments_final %>%
  filter(Microhabitat == "Mesic/Xeric")

light_mesic_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Light, 
                                  random = ~ animal + ID + n:ID, 
                                  family = "multinomial2", pedigree = meta_light_tree, 
                                  prior = priors, data = light_experiments_mesic_xeric,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

summary(light_mesic_xeric)
n_distinct(light_experiments_mesic_xeric$animal)

lambda <- light_mesic_xeric$VCV[,'animal']/
  (light_mesic_xeric$VCV[,'animal']+light_mesic_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Species with data on seed mass.

light_experiments_mass <- light_experiments_final %>%
  filter(!is.na(Dry_mass))

light_mass <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Light * scale(Dry_mass), 
                                        random = ~ animal + ID + n:ID, 
                                        family = "multinomial2", pedigree = meta_light_tree, 
                                        prior = priors, data = light_experiments_mass,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

summary(light_mass)
n_distinct(light_experiments_mass$animal)

lambda <- light_mass$VCV[,'animal']/
  (light_mass$VCV[,'animal']+light_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)
