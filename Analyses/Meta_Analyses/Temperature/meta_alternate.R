# The following code was designed to assess the effect of alternate temperatures on germination percentage and median germination time by performing meta-analyses through fitting binomial phylogenetic generalized mixed models with Bayesian estimation using Markov chain Monte Carlo (MCMC) method as implemented in the MCMCglmm R-package.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: April 17th, 2024.

#### 0. Loading packages (assuming they are already installed) ####

library(readr) # v. 2.1.5. for database loading.
library(dplyr) # v. 1.1.4 for database manipulation.
library(tidyr) # v. 1.3.0. (also) for database manipulation
library(tibble) # v. 3.2.1 (also) for database manipulation
library(forcats) # v. 1.0.0 for working with categorical variables.
library(picante) # v. 1.8.2 for loading phylogenetic trees
library(germinationmetrics) # v. 0.1.8.9000 to calculate t50
library(MCMCglmm) # v. 2.35 for generalized mixed models with Bayesian estimation using MCMC method.

#### 1. Setting working directory ####

setwd(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature"))

#### 2. Loading datasets ####

# For this script we will be using four data sets. First, we are loading a data frame that classifies the studies into the topics they studied. For this classification three columns are used: Factor, Major and Minor.

topics <- read_csv(str_c(git_find(),"/topics.csv"))

# Second, we are uploading a database with the germination data provided by the authors.

germination <- read_csv(str_c(git_find(),"/germination.csv"))

# Third, we are uploading the database of seed functional traits, which includes information on species distribution, microhabitats, growth form and seed mass.

traits <- read_csv(str_c(git_find(),"/traits.csv"))

# We are uploading a database with the species names standardized against the Leipzig Catalogue of Vascular Plants (Freiberg et al. 2020)

lcvp_names_checked <- read.csv(str_c(git_find(),"/lcvp_names_checked.csv"))

# We are also uploading the phylogenetic tree with all the species in our dataset.

full_tree <- read.tree(str_c(git_find(),"/Analyses/Phylogeny/outcrop_phylo.tre"))

#### 2. Getting the data from studies that assessed the effect of alternate temperatures ####

# The following code creates a data set of the germination experiments that tested the effect of light availability. The control treatment in all cases were seeds incubated under dark conditions.First we're getting the list of studies researching the effects of light availabilty on germination. Therefore, we are filtering the "topics" database looking for studies where "Major" is "Light" and "Minor" is "Availability".

alternate <- topics %>%
  filter(Minor == "Alternate")

# Then, we are getting the IDs of each of these studies.

alternate_studies <- unique(alternate$ID)
alternate_studies

# Now, we filter the germination database looking for:

alternate_experiments <- germination %>%
  filter(ID %in% alternate_studies, # Studies with whose IDs are in "light_studies",
         Photoperiod != "0",
         Temperature == "25" | grepl('/', Temperature),
         Chemical_Compound == "Control" | is.na(Chemical_Compound), # Experiments that only tested light effects.
         Storage_Time == "Control" | is.na(Storage_Time),
         Scarification == "Sandpaper" | is.na(Scarification),
         Gut_Passage == "Extracted seeds"| is.na(Gut_Passage),
         HeatShock_Temperature == "Control" | is.na(HeatShock_Temperature),
         Germination_Substrate == "Control" | is.na(Germination_Substrate),
         Smoke == "Control" | is.na(Smoke)) %>%
  dplyr::select(1:3, Temperature, GermSeeds, nSeeds) %>%
  left_join(traits) %>%
  dplyr::select(1:7, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits we're interested in.
  left_join(lcvp_names_checked) %>%
  filter(Growth_form == "Herb" | Growth_form == "Shrub",
         Temperature != "30/15") %>%
  mutate(Temperature = as.factor(Temperature))
  

levels(alternate_experiments$Temperature)

# to avoid the confounding effects from treatments with extremely low germination, we only included observation where either treatment (light and darkness) had >10% germination. For that purpose, we will create a new object with such observations

exclude_alternate_experiments <- alternate_experiments %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25/15" ~ "Alt25_15",
                                 Temperature == "30/20" ~ "Alt30_20",
                                 Temperature == "25" ~ "Control")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Alt30_20) & is.na(Alt25_15) | Control < 0.1 & Alt30_20 < 0.1 | Control < 0.1 & Alt25_15 < 0.1) %>%
  ungroup() %>%
  pivot_longer(cols = c(Control, Alt30_20, Alt25_15), names_to = "Temperature") %>%
  mutate(Temperature = case_when(Temperature == "Alt25_15" ~ "25/15",
                                 Temperature == "Alt30_20" ~ "30/20",
                                 Temperature == "Control" ~ "25")) %>%
  select(1:4)

alternate_experiments2 <- alternate_experiments %>%
  anti_join(exclude_alternate_experiments) %>%
  filter(Dormancy == "ND")

# Just to have an idea of the size of our dataset, we will look and how many species and studies we have.

n_distinct(alternate_experiments2$ID) # We have four studies
n_distinct(alternate_experiments2$Species_acceptedLCVP) # For a total of 13 species.

# Also, just to confirm which studies are being include in the meta-analysis, we will create the following object.

included_alternate_studies <- alternate_experiments2 %>% # Take the light_experiments2 dataset
  dplyr::select(ID) %>% # Select the 'ID" column
  distinct() # And keep only unique records (i.e., if the ID is repeated, keep only the first line).

included_alternate_studies

# Before proceeding to the next step, we need to input missing seed mass data. Thus, we'll be exporting the light_experiments2 dataset, to manually add the seed mass values.

write.csv(alternate_experiments2, file = "alternate_experiments_missing.csv")

# After imputing the missing seed mass data manually, we'll be loading the final dataset that we are going to use to fit the model. We'll be also adding a new column call n, that will allow us to control the effect of various observation between studies.

alternate_experiments_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/alternate_experiments_final.csv")) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()))

# Now, we will proceed to fit the models, first, by assessing the effect of 25/15 regimes.

alternate_experiments_25_15 <- alternate_experiments_final %>%
  filter(Condition != "30/20") %>%
  mutate(Condition = case_when(Condition == "25" ~ 0,
                               Condition == "25/15" ~ 1),
         Condition = as.factor(Condition))

meta_alternate25_15_names <- alternate_experiments_25_15 %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal") %>%
  names_for_pruning()

meta_alternate25_15_tree <- prune.sample(meta_alternate25_15_names, phylo = full_tree)
meta_alternate25_15_tree$node.label <- NULL

# Now, we will set the number of desired iterations based on De Villemereuil & Nakagawa (2014)

nite = 500000
nthi = 50
nbur = 50000

# Set priors based on https://doi.org/10.1017/S0960258517000332

priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

alternate25_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                random = ~ animal + ID + n:ID, 
                                family = "multinomial2", pedigree = meta_alternate25_15_tree, 
                                prior = priors, data = alternate_experiments_25_15,
                                nitt = nite, thin = nthi, burnin = nbur, 
                                verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                pr = FALSE, pl = FALSE)


summary(alternate25_all)

n_distinct(alternate_experiments_25_15$animal)

lambda <- alternate25_all$VCV[,'animal']/
  (alternate25_all$VCV[,'animal']+alternate25_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

alternate_experiments_25_15_mass <- alternate_experiments_25_15 %>%
  filter(!is.na(Dry_mass))

alternate25_mass <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition * scale(Dry_mass), 
                                 random = ~ animal + ID + n:ID, 
                                 family = "multinomial2", pedigree = meta_alternate25_15_tree, 
                                 prior = priors, data = alternate_experiments_25_15_mass,
                                 nitt = nite, thin = nthi, burnin = nbur, 
                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                 pr = FALSE, pl = FALSE)

summary(alternate25_mass)
n_distinct(alternate_experiments_25_15_mass$animal)

lambda <- alternate25_mass$VCV[,'animal']/
  (alternate25_mass$VCV[,'animal']+alternate25_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Now, we will proceed to fit the models, first, by assessing the effect of 30/20 regimes.

alternate_experiments_30_20 <- alternate_experiments_final %>%
  filter(Condition != "25/15") %>%
  mutate(Condition = case_when(Condition == "25" ~ 0,
                               Condition == "30/20" ~ 1),
         Condition = as.factor(Condition))

meta_alternate30_20_names <- alternate_experiments_30_20 %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal") %>%
  names_for_pruning()

meta_alternate30_20_tree <- prune.sample(meta_alternate30_20_names, phylo = full_tree)
meta_alternate30_20_tree$node.label <- NULL

alternate30_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                      random = ~ animal + ID + n:ID, 
                                      family = "multinomial2", pedigree = meta_alternate30_20_tree, 
                                      prior = priors, data = alternate_experiments_30_20,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)


summary(alternate30_all)

n_distinct(alternate_experiments_30_20$animal)

lambda <- alternate30_all$VCV[,'animal']/
  (alternate30_all$VCV[,'animal']+alternate30_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

alternate_experiments_30_20_herbs <- alternate_experiments_30_20 %>%
  filter(Growth_form == "Herb")

alternate30_herbs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                      random = ~ animal + ID + n:ID, 
                                      family = "multinomial2", pedigree = meta_alternate30_20_tree, 
                                      prior = priors, data = alternate_experiments_30_20_herbs,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

summary(alternate30_herbs)
n_distinct(alternate_experiments_30_20_herbs$animal)

lambda <- alternate30_herbs$VCV[,'animal']/
  (alternate30_herbs$VCV[,'animal']+alternate30_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

alternate_experiments_30_20_shrubs <- alternate_experiments_30_20 %>%
  filter(Growth_form == "Shrub")

alternate30_shrubs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "multinomial2", pedigree = meta_alternate30_20_tree, 
                                        prior = priors, data = alternate_experiments_30_20_shrubs,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)


summary(alternate30_shrubs)
n_distinct(alternate_experiments_30_20_shrubs$animal)

lambda <- alternate30_shrubs$VCV[,'animal']/
  (alternate30_shrubs$VCV[,'animal']+alternate30_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

alternate_experiments_30_20_mass <- alternate_experiments_30_20 %>%
  filter(!is.na(Dry_mass))

alternate30_mass <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition * scale(Dry_mass), 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_alternate30_20_tree, 
                                       prior = priors, data = alternate_experiments_30_20_mass,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)
beepr::beep()

summary(alternate30_mass)
n_distinct(alternate_experiments_30_20_mass$animal)

lambda <- alternate30_mass$VCV[,'animal']/
  (alternate30_mass$VCV[,'animal']+alternate30_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

## t50

alternate_experiments_time <- germination %>%
  filter(ID %in% alternate_studies, # Studies with whose IDs are in "light_studies",
         Photoperiod != "0",
         Temperature == "25" | grepl('/', Temperature),
         Chemical_Compound == "Control" | is.na(Chemical_Compound), # Experiments that only tested light effects.
         Storage_Time == "Control" | is.na(Storage_Time),
         Scarification == "Sandpaper" | is.na(Scarification),
         Gut_Passage == "Extracted seeds"| is.na(Gut_Passage),
         HeatShock_Temperature == "Control" | is.na(HeatShock_Temperature),
         Germination_Substrate == "Control" | is.na(Germination_Substrate),
         Smoke == "Control" | is.na(Smoke)) %>%
  dplyr::select(1:3, Temperature, GermSeeds, nSeeds, 26:55) %>%
  left_join(traits) %>%
  dplyr::select(1:36, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits we're interested in.
  left_join(lcvp_names_checked) %>%
  filter(Growth_form == "Herb" | Growth_form == "Shrub",
         Temperature != "30/15",
         !is.na(D1)) %>%
  mutate(Temperature = as.factor(Temperature))

exclude_alternate_experiments_time <- alternate_experiments_time %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25/15" ~ "Alt25_15",
                                 Temperature == "30/20" ~ "Alt30_20",
                                 Temperature == "25" ~ "Control")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Alt30_20) & is.na(Alt25_15) | !(Control > 0.1 & (Alt30_20 > 0.1 | Alt25_15 > 0.1))) %>%
  ungroup() %>%
  pivot_longer(cols = c(Control, Alt30_20, Alt25_15), names_to = "Temperature") %>%
  mutate(Temperature = case_when(Temperature == "Alt25_15" ~ "25/15",
                                 Temperature == "Alt30_20" ~ "30/20",
                                 Temperature == "Control" ~ "25")) %>%
  select(1:4)

alternate_experiments_time2 <- alternate_experiments_time %>%
  anti_join(exclude_alternate_experiments_time) %>%
  filter(Dormancy == "ND") %>%
  mutate(D0 = 0,
         n = as.factor(row_number())) %>%
  select(n, ID:nSeeds, D0, D1:Species_acceptedLCVP)

intervals30 <- paste0("D", 0:30)

alternate_experiments_t50 <- germination.indices(as.data.frame(alternate_experiments_time2),
                                                 total.seeds.col = "nSeeds",
                                                 counts.intervals.cols = intervals30,
                                                 intervals = 0:30,
                                                 t50 = TRUE,
                                                 EmergenceRateIndex = FALSE,
                                                 GermValue = FALSE) %>%
  dplyr::select(n:Species_acceptedLCVP, t50_Farooq)

### Saving the dataset to input seed mass data manually ####

write.csv(alternate_experiments_t50, file = "alternate_experiments_t50_missing.csv")

####

alternate_experiments_t50_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/alternate_experiments_t50_final.csv"))

alternate_experiments_t50_25 <- alternate_experiments_t50_final %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  filter(Condition != "30/20") %>%
  mutate(Condition = case_when(Condition == "25" ~ 0,
                               Condition == "25/15" ~ 1),
         Condition = as.factor(Condition))

alternate_experiments_t50_25_names <- alternate_experiments_t50_25 %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal") %>%
  names_for_pruning()

alternate_experiments_t50_25_tree <- prune.sample(alternate_experiments_t50_25_names, phylo = full_tree)
alternate_experiments_t50_25_tree$node.label <- NULL

alternate25_t50_all <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                      random = ~ animal + ID + n:ID, 
                                      family = "gaussian", pedigree = alternate_experiments_t50_25_tree, 
                                      prior = priors, data = alternate_experiments_t50_25,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()
summary(alternate25_t50_all)

n_distinct(alternate_experiments_t50_25$animal)

lambda <- alternate25_t50_all$VCV[,'animal']/
  (alternate25_t50_all$VCV[,'animal']+alternate25_t50_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

alternate_experiments_t50_25_mass <- alternate_experiments_t50_25 %>%
  filter(!is.na(Dry_mass))

alternate25_t50_mass <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition * scale(Dry_mass), 
                                          random = ~ animal + ID + n:ID, 
                                          family = "gaussian", pedigree = alternate_experiments_t50_25_tree, 
                                          prior = priors, data = alternate_experiments_t50_25_mass,
                                          nitt = nite, thin = nthi, burnin = nbur, 
                                          verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                          pr = FALSE, pl = FALSE)
beepr::beep()

summary(alternate25_t50_mass)

n_distinct(alternate_experiments_t50_25_mass$animal)

lambda <- alternate25_t50_mass$VCV[,'animal']/
  (alternate25_t50_mass$VCV[,'animal']+alternate25_t50_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

alternate_experiments_t50_30 <- alternate_experiments_t50_final %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  filter(Condition != "25/15") %>%
  mutate(Condition = case_when(Condition == "25" ~ 0,
                               Condition == "30/20" ~ 1),
         Condition = as.factor(Condition))

alternate_experiments_t50_30_names <- alternate_experiments_t50_30 %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal") %>%
  names_for_pruning()

alternate_experiments_t50_30_tree <- prune.sample(alternate_experiments_t50_30_names, phylo = full_tree)
alternate_experiments_t50_30_tree$node.label <- NULL

alternate30_t50_all <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                          random = ~ animal + ID + n:ID, 
                                          family = "gaussian", pedigree = alternate_experiments_t50_30_tree, 
                                          prior = priors, data = alternate_experiments_t50_30,
                                          nitt = nite, thin = nthi, burnin = nbur, 
                                          verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                          pr = FALSE, pl = FALSE)

beepr::beep()

summary(alternate30_t50_all)

n_distinct(alternate_experiments_t50_30$animal)

lambda <- alternate30_t50_all$VCV[,'animal']/
  (alternate30_t50_all$VCV[,'animal']+alternate30_t50_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

alternate_experiments_t50_30_herbs <- alternate_experiments_t50_30 %>%
  filter(Growth_form == "Herb")

alternate30_t50_herbs <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                          random = ~ animal + ID + n:ID, 
                                          family = "gaussian", pedigree = alternate_experiments_t50_30_tree, 
                                          prior = priors, data = alternate_experiments_t50_30_herbs,
                                          nitt = nite, thin = nthi, burnin = nbur, 
                                          verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                          pr = FALSE, pl = FALSE)

beepr::beep()

summary(alternate30_t50_herbs)

n_distinct(alternate_experiments_t50_30_herbs$animal)

lambda <- alternate30_t50_herbs$VCV[,'animal']/
  (alternate30_t50_herbs$VCV[,'animal']+alternate30_t50_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

alternate_experiments_t50_30_shrubs <- alternate_experiments_t50_30 %>%
  filter(Growth_form == "Shrub")

alternate30_t50_shrubs <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = alternate_experiments_t50_30_tree, 
                                            prior = priors, data = alternate_experiments_t50_30_shrubs,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(alternate30_t50_shrubs)

n_distinct(alternate_experiments_t50_30_shrubs$animal)

lambda <- alternate30_t50_shrubs$VCV[,'animal']/
  (alternate30_t50_shrubs$VCV[,'animal']+alternate30_t50_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

alternate_experiments_t50_30_mass <- alternate_experiments_t50_30 %>%
  filter(!is.na(Dry_mass))

alternate30_t50_mass <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "gaussian", pedigree = alternate_experiments_t50_30_tree, 
                                             prior = priors, data = alternate_experiments_t50_30_mass,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(alternate30_t50_mass)

n_distinct(alternate_experiments_t50_30_mass$animal)

lambda <- alternate30_t50_mass$VCV[,'animal']/
  (alternate30_t50_mass$VCV[,'animal']+alternate30_t50_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)
