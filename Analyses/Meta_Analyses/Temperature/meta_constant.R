# The following code was designed to assess the effect of constant temperatures on germination percentage and median germination time by performing meta-analyses through fitting binomial phylogenetic generalized mixed models with Bayesian estimation using Markov chain Monte Carlo (MCMC) method as implemented in the MCMCglmm R-package.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: April 25th, 2024.

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

#### 3. Getting the data from studies that assessed the effect of temperature at 10 °C ####

# The following code looks to create a data set of the germination experiments that tested the effect of constant temperatures. The control treatment in all cases were seeds incubated under 25 ºC.

# First we're getting the list of studies researching the effects of constant temperatures on germination. Therefore, we are filtering the "topics" database looking for studies where "Major" is "Temperature" and "Minor" is "Constant".

constant <- topics %>%
  filter(Major == "Temperature",
         Minor == "Constant")

# Then, we are getting the IDs of each of these studies.

constant_studies <- unique(constant$ID)
constant_studies

# Now, we're filtering the germination database looking for:

constant_experiments <- germination %>%
  filter(ID %in% constant_studies, # Studies with whose IDs are in "constant_studies",
         Photoperiod != 0,
         !grepl('/', Temperature),
         Chemical_Compound == "Control" | is.na(Chemical_Compound),
         Storage_Time == "Control" | is.na(Storage_Time),
         Scarification == "Sandpaper" | is.na(Scarification),
         Germination_Substrate == "Control"| is.na(Germination_Substrate),
         HeatShock_Temperature == "Control" | is.na(HeatShock_Temperature),
         Smoke == "Control" | is.na(Smoke),
         Hypoxia == "Control" | is.na(Hypoxia)) %>%
  select(1:3, Temperature, GermSeeds, nSeeds) %>%
  left_join(traits) %>%
  dplyr::select(1:6, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits we're interested in.
  left_join(lcvp_names_checked) %>%
  filter(Growth_form == "Herb" | Growth_form == "Shrub",
         !(Dormancy == "D" & !(ID %in% c(14, 78))), 
                Dormancy != "NC") %>%
  mutate(Temperature = as.factor(Temperature))

constant_experiments10 <- constant_experiments %>%
  filter(Temperature == "10" | Temperature == "25")

exclude_constant_experiments10 <- constant_experiments %>%
  filter(Temperature == "10" | Temperature == "25") %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25" ~ "Control",
                                 Temperature == "10" ~ "Treatment")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 & Treatment < 0.1)) %>%
  select(1:3)

constant_experiments10_2 <- constant_experiments10 %>%
  anti_join(exclude_constant_experiments10)
         
# Creating the dataset for adding missing seed mass values

write.csv(constant_experiments10_2, file = "constant_experiments10_missing.csv")

# Loading the dataset with the inputed mass data.

constant_experiments10_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/constant_experiments10_final.csv")) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Condition = case_when(Condition == "10" ~ 1,
                               Condition == "25" ~ 0))

meta_constant10_names <- constant_experiments10_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_constant10_names <- names_for_pruning(meta_constant10_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_constant10_tree <- prune.sample(meta_constant10_names, phylo = full_tree)
meta_constant10$node.label <- NULL

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

constant10_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                random = ~ animal + ID + n:ID, 
                                family = "multinomial2", pedigree = meta_constant10_tree, 
                                prior = priors, data = constant_experiments10_final,
                                nitt = nite, thin = nthi, burnin = nbur, 
                                verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant10_all)
plot(constant10_all)

n_distinct(constant_experiments10_final$animal)

lambda <- constant10_all$VCV[,'animal']/
  (constant10_all$VCV[,'animal']+constant10_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments10_restricted <- constant_experiments10_final %>%
  filter(Distribution == "Restricted")

constant10_restricted <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                     random = ~ animal + ID + n:ID, 
                                     family = "multinomial2", pedigree = meta_constant10_tree, 
                                     prior = priors, data = constant_experiments10_restricted,
                                     nitt = nite, thin = nthi, burnin = nbur, 
                                     verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                     pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant10_restricted)
n_distinct(constant_experiments10_restricted$animal)

lambda <- constant10_restricted$VCV[,'animal']/
  (constant10_restricted$VCV[,'animal']+constant10_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments10_widespread <- constant_experiments10_final %>%
  filter(Distribution == "Widespread")

priors2 <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

constant10_widespread <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + n, # Here we don't control for ID, since we only have                                               one study 
                                            family = "multinomial2", pedigree = meta_constant10_tree, 
                                            prior = priors2, data = constant_experiments10_widespread,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant10_widespread)
n_distinct(constant_experiments10_widespread$animal)

lambda <- constant10_widespread$VCV[,'animal']/
  (constant10_widespread$VCV[,'animal']+constant10_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments10_mass <- constant_experiments10_final %>%
  filter(!is.na(Dry_mass))

constant10_mass <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition * scale(Dry_mass), 
                                            random = ~ animal + ID + n:ID, 
                                            family = "multinomial2", pedigree = meta_constant10_tree, 
                                            prior = priors, data = constant_experiments10_mass,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant10_mass)

n_distinct(constant_experiments10_mass$animal)

lambda <- constant10_mass$VCV[,'animal']/
  (constant10_mass$VCV[,'animal']+constant10_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

#### 4. 15 °C ####

constant_experiments15 <- constant_experiments %>%
  filter(Temperature == "15" | Temperature == "25")

exclude_constant_experiments15 <- constant_experiments %>%
  filter(Temperature == "15" | Temperature == "25") %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25" ~ "Control",
                                 Temperature == "15" ~ "Treatment")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 & Treatment < 0.1)) %>%
  select(1:3)

constant_experiments15_2 <- constant_experiments15 %>%
  anti_join(exclude_constant_experiments15)

# Creating the dataset for adding missing seed mass values

write.csv(constant_experiments15_2, file = "constant_experiments15_missing.csv")

# Loading the dataset with the inputed mass data.

constant_experiments15_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/constant_experiments15_final.csv")) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Condition = case_when(Condition == "15" ~ 1,
                               Condition == "25" ~ 0),
         Condition = as.factor(Condition))

meta_constant15_names <- constant_experiments15_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_constant15_names <- names_for_pruning(meta_constant15_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_constant15_tree <- prune.sample(meta_constant15_names, phylo = full_tree)
meta_constant15_tree$node.label <- NULL

constant15_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                     random = ~ animal + ID + n:ID, 
                                     family = "multinomial2", pedigree = meta_constant15_tree, 
                                     prior = priors, data = constant_experiments15_final,
                                     nitt = nite, thin = nthi, burnin = nbur, 
                                     verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                     pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant15_all)
n_distinct(constant_experiments15_final$animal)
nrow(constant_experiments15_final)

lambda <- constant15_all$VCV[,'animal']/
  (constant15_all$VCV[,'animal']+constant15_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments15_herbs <- constant_experiments15_final %>%
  filter(Growth_form == "Herb")

constant15_herbs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                     random = ~ animal + ID + n:ID, 
                                     family = "multinomial2", pedigree = meta_constant15_tree, 
                                     prior = priors, data = constant_experiments15_herbs,
                                     nitt = nite, thin = nthi, burnin = nbur, 
                                     verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                     pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant15_herbs)
n_distinct(constant_experiments15_herbs$animal)
nrow(constant_experiments15_herbs)

lambda <- constant15_herbs$VCV[,'animal']/
  (constant15_herbs$VCV[,'animal']+constant15_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments15_shrubs <- constant_experiments15_final %>%
  filter(Growth_form == "Shrub")

constant15_shrubs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant15_tree, 
                                       prior = priors, data = constant_experiments15_shrubs,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant15_shrubs)
n_distinct(constant_experiments15_shrubs$animal)
nrow(constant_experiments15_shrubs)

lambda <- constant15_shrubs$VCV[,'animal']/
  (constant15_shrubs$VCV[,'animal']+constant15_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments15_restricted <- constant_experiments15_final %>%
  filter(Distribution == "Restricted")

n_distinct(constant_experiments15_shrubs$animal)
nrow(constant_experiments15_shrubs)

constant15_restricted <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "multinomial2", pedigree = meta_constant15_tree, 
                                        prior = priors, data = constant_experiments15_restricted,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant15_restricted)


lambda <- constant15_restricted$VCV[,'animal']/
  (constant15_restricted$VCV[,'animal']+constant15_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments15_widespread <- constant_experiments15_final %>%
  filter(Distribution == "Widespread")

constant15_widespread <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "multinomial2", pedigree = meta_constant15_tree, 
                                            prior = priors, data = constant_experiments15_widespread,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant15_widespread)
n_distinct(constant_experiments15_widespread$animal)
nrow(constant_experiments15_widespread)

lambda <- constant15_widespread$VCV[,'animal']/
  (constant15_widespread$VCV[,'animal']+constant15_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments15_mesic <- constant_experiments15_final %>%
  filter(Microhabitat == "Mesic")

constant15_mesic <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "multinomial2", pedigree = meta_constant15_tree, 
                                            prior = priors, data = constant_experiments15_mesic,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant15_mesic)
n_distinct(constant_experiments15_mesic$animal)
nrow(constant_experiments15_mesic)

lambda <- constant15_mesic$VCV[,'animal']/
  (constant15_mesic$VCV[,'animal']+constant15_mesic$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments15_mesic_xeric <- constant_experiments15_final %>%
  filter(Microhabitat == "Mesic/Xeric")

constant15_mesic_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant15_tree, 
                                       prior = priors, data = constant_experiments15_mesic_xeric,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant15_mesic_xeric)
n_distinct(constant_experiments15_mesic_xeric$animal)
nrow(constant_experiments15_mesic_xeric)

lambda <- constant15_mesic_xeric$VCV[,'animal']/
  (constant15_mesic_xeric$VCV[,'animal']+constant15_mesic_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments15_xeric <- constant_experiments15_final %>%
  filter(Microhabitat == "Xeric")

constant15_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "multinomial2", pedigree = meta_constant15_tree, 
                                             prior = priors, data = constant_experiments15_xeric,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant15_xeric)
n_distinct(constant_experiments15_xeric$animal)
nrow(constant_experiments15_xeric)

lambda <- constant15_xeric$VCV[,'animal']/
  (constant15_xeric$VCV[,'animal']+constant15_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments15_mass <- constant_experiments15_final %>%
  filter(!is.na(Dry_mass))

constant15_mass <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition * scale(Dry_mass), 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant15_tree, 
                                       prior = priors, data = constant_experiments15_mass,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant15_mass)
n_distinct(constant_experiments15_mass$animal)
nrow(constant_experiments15_mass)

lambda <- constant15_mass$VCV[,'animal']/
  (constant15_mass$VCV[,'animal']+constant15_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

#t50

constant_experiments_time <- germination %>%
  filter(ID %in% constant_studies, # Studies with whose IDs are in "constant_studies",
         Photoperiod != 0,
         !grepl('/', Temperature),
         Chemical_Compound == "Control" | is.na(Chemical_Compound),
         Storage_Time == "Control" | is.na(Storage_Time),
         Scarification == "Sandpaper" | is.na(Scarification),
         Germination_Substrate == "Control"| is.na(Germination_Substrate),
         HeatShock_Temperature == "Control" | is.na(HeatShock_Temperature),
         Smoke == "Control" | is.na(Smoke),
         Hypoxia == "Control" | is.na(Hypoxia)) %>%
  select(1:3, Temperature, GermSeeds, nSeeds, 26:55) %>%
  left_join(traits) %>%
  dplyr::select(1:36, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits we're interested in.
  left_join(lcvp_names_checked) %>%
  filter(!is.na(D1),
         Growth_form == "Herb" | Growth_form == "Shrub",
         !(Dormancy == "D" & !(ID %in% c(14, 78))), 
         Dormancy != "NC") %>%
  mutate(Temperature = as.factor(Temperature))

constant_experiments_time15 <- constant_experiments_time %>%
  filter(Temperature == "15" | Temperature == "25")

exclude_constant_experiments_time15 <- constant_experiments_time %>%
  filter(Temperature == "15" | Temperature == "25") %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25" ~ "Control",
                                 Temperature == "15" ~ "Treatment")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 | Treatment < 0.1)) %>%
  select(1:3)

constant_experiments_time15_2 <- constant_experiments_time15 %>%
  anti_join(exclude_constant_experiments_time15) %>%
  mutate(D0 = 0,
         n = as.factor(row_number())) %>%
  select(n, ID:nSeeds, D0, D1:Species_acceptedLCVP)

constant_experiments_time15_t50 <- germination.indices(as.data.frame(constant_experiments_time15_2),
                                                       total.seeds.col = "nSeeds",
                                                       counts.intervals.cols = intervals30,
                                                       intervals = 1:31,
                                                       t50 = TRUE,
                                                       EmergenceRateIndex = FALSE,
                                                       PeakGermPercent = FALSE,
                                                       PeakGermTime = FALSE,
                                                       GermValue = FALSE) %>%
  dplyr::select(n:Temperature, Growth_form:Species_acceptedLCVP, t50_Farooq)

write.csv(constant_experiments_time15_t50, file = "constant_experiments_time15_missing.csv")

constant_experiments_time15_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/constant_experiments_time15_final.csv")) %>%
  filter(!is.na(Species_acceptedLCVP)) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Condition = case_when(Condition == "15" ~ 1,
                               Condition == "25" ~ 0),
         Condition = as.factor(Condition))

meta_constant_time15_names <- constant_experiments_time15_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_constant_time15_names <- names_for_pruning(meta_constant_time15_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_constant_time15_tree <- prune.sample(meta_constant_time15_names, phylo = full_tree)
meta_constant_time15_tree$node.label <- NULL

constant_time15_all <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                          random = ~ animal + ID + n:ID, 
                                          family = "gaussian", pedigree = meta_constant_time15_tree, 
                                          prior = priors, data = constant_experiments_time15_final,
                                          nitt = nite, thin = nthi, burnin = nbur, 
                                          verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                          pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time15_all)
n_distinct(constant_experiments_time15_final$animal)
nrow(constant_experiments_time15_final)

lambda <- constant_time15_all$VCV[,'animal']/
  (constant_time15_all$VCV[,'animal']+constant_time15_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time15_herbs <- constant_experiments_time15_final %>%
  filter(Growth_form == "Herb")

constant_time15_herbs <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time15_tree, 
                                            prior = priors, data = constant_experiments_time15_herbs,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time15_herbs)
n_distinct(constant_experiments_time15_herbs$animal)
nrow(constant_experiments_time15_herbs)

lambda <- constant_time20_herbs$VCV[,'animal']/
  (constant_time20_herbs$VCV[,'animal']+constant_time20_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time15_shrubs <- constant_experiments_time15_final %>%
  filter(Growth_form == "Shrub")

constant_time15_shrubs <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "gaussian", pedigree = meta_constant_time15_tree, 
                                             prior = priors, data = constant_experiments_time15_shrubs,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time15_shrubs)
n_distinct(constant_experiments_time15_shrubs$animal)
nrow(constant_experiments_time15_shrubs)

lambda <- constant_time15_shrubs$VCV[,'animal']/
  (constant_time15_shrubs$VCV[,'animal']+constant_time15_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time15_restricted <- constant_experiments_time15_final %>%
  filter(Distribution == "Restricted")

constant_time15_restricted <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                 random = ~ animal + ID + n:ID, 
                                                 family = "gaussian", pedigree = meta_constant_time15_tree, 
                                                 prior = priors, data = constant_experiments_time15_restricted,
                                                 nitt = nite, thin = nthi, burnin = nbur, 
                                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                 pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time15_restricted)
n_distinct(constant_experiments_time15_restricted$animal)
nrow(constant_experiments_time15_restricted)

lambda <- constant_time15_restricted$VCV[,'animal']/
  (constant_time15_restricted$VCV[,'animal']+constant_time15_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time15_widespread <- constant_experiments_time15_final %>%
  filter(Distribution == "Widespread")

constant_time15_widespread <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                 random = ~ animal + ID + n:ID, 
                                                 family = "gaussian", pedigree = meta_constant_time15_tree, 
                                                 prior = priors, data = constant_experiments_time15_widespread,
                                                 nitt = nite, thin = nthi, burnin = nbur, 
                                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                 pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time15_widespread)
n_distinct(constant_experiments_time15_widespread$animal)
nrow(constant_experiments_time15_widespread)

lambda <- constant_time15_widespread$VCV[,'animal']/
  (constant_time15_widespread$VCV[,'animal']+constant_time15_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time15_mesic <- constant_experiments_time15_final %>%
  filter(Microhabitat == "Mesic")

constant_time15_mesic <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time15_tree,
                                            prior = priors, data = constant_experiments_time15_mesic,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time15_mesic)
n_distinct(constant_experiments_time15_mesic$animal)
nrow(constant_experiments_time15_mesic)

lambda <- constant_time15_mesic$VCV[,'animal']/
  (constant_time15_mesic$VCV[,'animal']+constant_time15_mesic$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time15_mesic_xeric <- constant_experiments_time15_final %>%
  filter(Microhabitat == "Mesic/Xeric")

constant_time15_mesic_xeric <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                  random = ~ animal + ID + n:ID, 
                                                  family = "gaussian", pedigree = meta_constant_time15_tree,
                                                  prior = priors, data = constant_experiments_time15_mesic_xeric,
                                                  nitt = nite, thin = nthi, burnin = nbur, 
                                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time15_mesic_xeric)
n_distinct(constant_experiments_time15_mesic_xeric$animal)
nrow(constant_experiments_time15_mesic_xeric)


lambda <- constant_time15_mesic_xeric$VCV[,'animal']/
  (constant_time15_mesic_xeric$VCV[,'animal']+constant_time15_mesic_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time15_xeric <- constant_experiments_time15_final %>%
  filter(Microhabitat == "Xeric")

constant_time15_xeric <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time15_tree, 
                                            prior = priors, data = constant_experiments_time15_xeric,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time15_xeric)
n_distinct(constant_experiments15_xeric$animal)
nrow(constant_experiments_time15_xeric)

lambda <- constant_time15_xeric$VCV[,'animal']/
  (constant_time15_xeric$VCV[,'animal']+constant_time15_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time15_mass <- constant_experiments_time15_final %>%
  filter(!is.na(Dry_mass))

constant_time15_mass <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition * scale(Dry_mass), 
                                           random = ~ animal + ID + n:ID, 
                                           family = "gaussian", pedigree = meta_constant_time15_tree, 
                                           prior = priors, data = constant_experiments_time15_mass,
                                           nitt = nite, thin = nthi, burnin = nbur, 
                                           verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                           pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time15_mass)
n_distinct(constant_experiments_time15_mass$animal)
nrow(constant_experiments_time15_mass)

lambda <- constant_time15_mass$VCV[,'animal']/
  (constant_time15_mass$VCV[,'animal']+constant_time15_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

#### 5. 20 °C ####

constant_experiments20 <- constant_experiments %>%
  filter(Temperature == "20" | Temperature == "25")

exclude_constant_experiments20 <- constant_experiments %>%
  filter(Temperature == "20" | Temperature == "25") %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25" ~ "Control",
                                 Temperature == "20" ~ "Treatment")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 & Treatment < 0.1)) %>%
  select(1:3)

constant_experiments20_2 <- constant_experiments20 %>%
  anti_join(exclude_constant_experiments20)

# Creating the dataset for adding missing seed mass values

write.csv(constant_experiments20_2, file = "constant_experiments20_missing.csv")

# Loading the dataset with the inputed mass data.

constant_experiments20_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/constant_experiments20_final.csv")) %>%
  filter(!is.na(Species_acceptedLCVP)) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Condition = case_when(Condition == "20" ~ 1,
                               Condition == "25" ~ 0),
         Condition = as.factor(Condition))

meta_constant20_names <- constant_experiments20_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_constant20_names <- names_for_pruning(meta_constant20_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_constant20_tree <- prune.sample(meta_constant20_names, phylo = full_tree)
meta_constant20_tree$node.label <- NULL

constant20_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                     random = ~ animal + ID + n:ID, 
                                     family = "multinomial2", pedigree = meta_constant20_tree, 
                                     prior = priors, data = constant_experiments20_final,
                                     nitt = nite, thin = nthi, burnin = nbur, 
                                     verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                     pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant20_all)
n_distinct(constant_experiments20_final$animal)
nrow(constant_experiments20_final)


lambda <- constant20_all$VCV[,'animal']/
  (constant20_all$VCV[,'animal']+constant20_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments20_herbs <- constant_experiments20_final %>%
  filter(Growth_form == "Herb")

constant20_herbs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant20_tree, 
                                       prior = priors, data = constant_experiments20_herbs,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant20_herbs)
n_distinct(constant_experiments20_herbs$animal)
nrow(constant_experiments20_herbs)

lambda <- constant20_herbs$VCV[,'animal']/
  (constant20_herbs$VCV[,'animal']+constant20_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments20_shrubs <- constant_experiments20_final %>%
  filter(Growth_form == "Shrub")

constant20_shrubs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "multinomial2", pedigree = meta_constant20_tree, 
                                        prior = priors, data = constant_experiments20_shrubs,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant20_shrubs)
n_distinct(constant_experiments20_shrubs$animal)
nrow(constant_experiments20_shrubs)

lambda <- constant20_shrubs$VCV[,'animal']/
  (constant20_shrubs$VCV[,'animal']+constant20_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments20_restricted <- constant_experiments20_final %>%
  filter(Distribution == "Restricted")

constant20_restricted <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "multinomial2", pedigree = meta_constant20_tree, 
                                            prior = priors, data = constant_experiments20_restricted,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant20_restricted)
n_distinct(constant_experiments20_restricted$animal)
nrow(constant_experiments20_restricted)


lambda <- constant20_restricted$VCV[,'animal']/
  (constant20_restricted$VCV[,'animal']+constant20_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments20_widespread <- constant_experiments20_final %>%
  filter(Distribution == "Widespread")

constant20_widespread <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "multinomial2", pedigree = meta_constant20_tree, 
                                            prior = priors, data = constant_experiments20_widespread,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant20_widespread)
n_distinct(constant_experiments20_widespread$animal)
nrow(constant_experiments20_widespread)

lambda <- constant20_widespread$VCV[,'animal']/
  (constant20_widespread$VCV[,'animal']+constant20_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments20_mesic <- constant_experiments20_final %>%
  filter(Microhabitat == "Mesic")

constant20_mesic <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant20_tree, 
                                       prior = priors, data = constant_experiments20_mesic,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant20_mesic)
n_distinct(constant_experiments20_mesic$animal)
nrow(constant_experiments20_mesic)

lambda <- constant20_mesic$VCV[,'animal']/
  (constant20_mesic$VCV[,'animal']+constant20_mesic$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments20_mesic_xeric <- constant_experiments20_final %>%
  filter(Microhabitat == "Mesic/Xeric")

constant20_mesic_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "multinomial2", pedigree = meta_constant20_tree, 
                                             prior = priors, data = constant_experiments20_mesic_xeric,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant20_mesic_xeric)
n_distinct(constant_experiments20_mesic_xeric$animal)
nrow(constant_experiments20_mesic_xeric)

lambda <- constant20_mesic_xeric$VCV[,'animal']/
  (constant20_mesic_xeric$VCV[,'animal']+constant20_mesic_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments20_xeric <- constant_experiments20_final %>%
  filter(Microhabitat == "Xeric")

constant20_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant20_tree, 
                                       prior = priors, data = constant_experiments20_xeric,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant20_xeric)
n_distinct(constant_experiments20_xeric$animal)
nrow(constant_experiments20_xeric)

lambda <- constant20_xeric$VCV[,'animal']/
  (constant20_xeric$VCV[,'animal']+constant20_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments20_mass <- constant_experiments20_final %>%
  filter(!is.na(Dry_mass))

constant20_mass <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition * scale(Dry_mass), 
                                      random = ~ animal + ID + n:ID, 
                                      family = "multinomial2", pedigree = meta_constant20_tree, 
                                      prior = priors, data = constant_experiments20_mass,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant20_mass)
n_distinct(constant_experiments20_xeric$animal)
nrow(constant_experiments20_mass)

lambda <- constant20_mass$VCV[,'animal']/
  (constant20_mass$VCV[,'animal']+constant20_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

#t50

constant_experiments_time20 <- constant_experiments_time %>%
  filter(Temperature == "20" | Temperature == "25")

exclude_constant_experiments_time20 <- constant_experiments_time %>%
  filter(Temperature == "20" | Temperature == "25") %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25" ~ "Control",
                                 Temperature == "20" ~ "Treatment")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 | Treatment < 0.1)) %>%
  select(1:3)

constant_experiments_time20_2 <- constant_experiments_time20 %>%
  anti_join(exclude_constant_experiments_time20) %>%
  mutate(D0 = 0,
         n = as.factor(row_number())) %>%
  select(n, ID:nSeeds, D0, D1:Species_acceptedLCVP)

length(intervals30)

constant_experiments_time20_t50 <- germination.indices(as.data.frame(constant_experiments_time20_2),
                                                 total.seeds.col = "nSeeds",
                                                 counts.intervals.cols = intervals30,
                                                 intervals = 1:31,
                                                 t50 = TRUE,
                                                 EmergenceRateIndex = FALSE,
                                                 PeakGermPercent = FALSE,
                                                 PeakGermTime = FALSE,
                                                 GermValue = FALSE) %>%
  dplyr::select(n:Temperature, Growth_form:Species_acceptedLCVP, t50_Farooq)

write.csv(constant_experiments_time20_t50, file = "constant_experiments_time20_missing.csv")

constant_experiments_time20_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/constant_experiments_time20_final.csv")) %>%
  filter(!is.na(Species_acceptedLCVP)) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Condition = case_when(Condition == "20" ~ 1,
                               Condition == "25" ~ 0),
         Condition = as.factor(Condition))

meta_constant_time20_names <- constant_experiments_time20_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_constant_time20_names <- names_for_pruning(meta_constant_time20_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_constant_time20_tree <- prune.sample(meta_constant_time20_names, phylo = full_tree)
meta_constant_time20_tree$node.label <- NULL

constant_time20_all <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                          random = ~ animal + ID + n:ID, 
                                          family = "gaussian", pedigree = meta_constant_time20_tree, 
                                          prior = priors, data = constant_experiments_time20_final,
                                          nitt = nite, thin = nthi, burnin = nbur, 
                                          verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                          pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time20_all)
n_distinct(constant_experiments_time20_final$animal)
nrow(constant_experiments_time20_final)

lambda <- constant_time20_all$VCV[,'animal']/
  (constant_time20_all$VCV[,'animal']+constant_time20_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time20_herbs <- constant_experiments_time20_final %>%
  filter(Growth_form == "Herb")

constant_time20_herbs <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "gaussian", pedigree = meta_constant_time20_tree, 
                                       prior = priors, data = constant_experiments_time20_herbs,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time20_herbs)
n_distinct(constant_experiments_time20_herbs$animal)
nrow(constant_experiments_time20_herbs)

lambda <- constant_time20_herbs$VCV[,'animal']/
  (constant_time20_herbs$VCV[,'animal']+constant_time20_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time20_shrubs <- constant_experiments_time20_final %>%
  filter(Growth_form == "Shrub")

constant_time20_shrubs <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "gaussian", pedigree = meta_constant_time20_tree, 
                                             prior = priors, data = constant_experiments_time20_shrubs,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time20_shrubs)
n_distinct(constant_experiments_time20_shrubs$animal)
nrow(constant_experiments_time20_shrubs)

lambda <- constant_time20_shrubs$VCV[,'animal']/
  (constant_time20_shrubs$VCV[,'animal']+constant_time20_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time20_restricted <- constant_experiments_time20_final %>%
  filter(Distribution == "Restricted")

constant_time20_restricted <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                 random = ~ animal + ID + n:ID, 
                                                 family = "gaussian", pedigree = meta_constant_time20_tree, 
                                                 prior = priors, data = constant_experiments_time20_restricted,
                                                 nitt = nite, thin = nthi, burnin = nbur, 
                                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                 pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time20_restricted)
n_distinct(constant_experiments_time20_restricted$animal)
nrow(constant_experiments_time20_restricted)

lambda <- constant_time20_restricted$VCV[,'animal']/
  (constant_time20_restricted$VCV[,'animal']+constant_time20_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time20_widespread <- constant_experiments_time20_final %>%
  filter(Distribution == "Widespread")

constant_time20_widespread <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                 random = ~ animal + ID + n:ID, 
                                                 family = "gaussian", pedigree = meta_constant_time20_tree, 
                                                 prior = priors, data = constant_experiments_time20_widespread,
                                                 nitt = nite, thin = nthi, burnin = nbur, 
                                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                 pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time20_widespread)
n_distinct(constant_experiments_time20_widespread$animal)
nrow(constant_experiments_time20_widespread)

lambda <- constant_time20_widespread$VCV[,'animal']/
  (constant_time20_widespread$VCV[,'animal']+constant_time20_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time20_mesic <- constant_experiments_time20_final %>%
  filter(Microhabitat == "Mesic")

constant_time20_mesic <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time20_tree,
                                            prior = priors, data = constant_experiments_time20_mesic,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time20_mesic)
n_distinct(constant_experiments_time20_mesic$animal)
nrow(constant_experiments_time20_mesic)

lambda <- constant_time20_mesic$VCV[,'animal']/
  (constant_time20_mesic$VCV[,'animal']+constant_time20_mesic$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time20_mesic_xeric <- constant_experiments_time20_final %>%
  filter(Microhabitat == "Mesic/Xeric")

constant_time20_mesic_xeric <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                  random = ~ animal + ID + n:ID, 
                                                  family = "gaussian", pedigree = meta_constant_time20_tree,
                                                  prior = priors, data = constant_experiments_time20_mesic_xeric,
                                                  nitt = nite, thin = nthi, burnin = nbur, 
                                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time20_mesic_xeric)
n_distinct(constant_experiments_time20_mesic_xeric$animal)
nrow(constant_experiments_time20_mesic_xeric)

lambda <- constant_time20_mesic_xeric$VCV[,'animal']/
  (constant_time20_mesic_xeric$VCV[,'animal']+constant_time20_mesic_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time20_xeric <- constant_experiments_time20_final %>%
  filter(Microhabitat == "Xeric")

constant_time20_xeric <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time20_tree, 
                                            prior = priors, data = constant_experiments_time20_xeric,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time20_xeric)
n_distinct(constant_experiments20_xeric$animal)
nrow(constant_experiments_time20_xeric)

lambda <- constant_time20_xeric$VCV[,'animal']/
  (constant_time20_xeric$VCV[,'animal']+constant_time20_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time20_mass <- constant_experiments_time20_final %>%
  filter(!is.na(Dry_mass))

constant_time20_mass <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition * scale(Dry_mass), 
                                      random = ~ animal + ID + n:ID, 
                                      family = "gaussian", pedigree = meta_constant_time20_tree, 
                                      prior = priors, data = constant_experiments_time20_mass,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time20_mass)
n_distinct(constant_experiments_time20_mass$animal)
nrow(constant_experiments_time20_mass)

lambda <- constant_time20_mass$VCV[,'animal']/
  (constant_time20_mass$VCV[,'animal']+constant_time20_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)


#### 6. 30 °C ####

constant_experiments30 <- constant_experiments %>%
  filter(Temperature == "30" | Temperature == "25")

exclude_constant_experiments30 <- constant_experiments %>%
  filter(Temperature == "30" | Temperature == "25") %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25" ~ "Control",
                                 Temperature == "30" ~ "Treatment")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 & Treatment < 0.1)) %>%
  select(1:3)

constant_experiments30_2 <- constant_experiments30 %>%
  anti_join(exclude_constant_experiments30)

# Creating the dataset for adding missing seed mass values

write.csv(constant_experiments30_2, file = "constant_experiments30_missing.csv")

# Loading the dataset with the inputed mass data.

constant_experiments30_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/constant_experiments30_final.csv")) %>%
  filter(!is.na(Species_acceptedLCVP)) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Condition = case_when(Condition == "30" ~ 1,
                               Condition == "25" ~ 0),
         Condition = as.factor(Condition))

meta_constant30_names <- constant_experiments30_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_constant30_names <- names_for_pruning(meta_constant30_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_constant30_tree <- prune.sample(meta_constant30_names, phylo = full_tree)
meta_constant30_tree$node.label <- NULL

constant30_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                     random = ~ animal + ID + n:ID, 
                                     family = "multinomial2", pedigree = meta_constant30_tree, 
                                     prior = priors, data = constant_experiments30_final,
                                     nitt = nite, thin = nthi, burnin = nbur, 
                                     verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                     pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant30_all)
n_distinct(constant_experiments30_final$animal)
nrow(constant_experiments30_final)

lambda <- constant30_all$VCV[,'animal']/
  (constant30_all$VCV[,'animal']+constant30_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments30_herbs <- constant_experiments30_final %>%
  filter(Growth_form == "Herb")

constant30_herbs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant30_tree, 
                                       prior = priors, data = constant_experiments30_herbs,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant30_herbs)
n_distinct(constant_experiments30_herbs$animal)
nrow(constant_experiments30_herbs)

lambda <- constant30_herbs$VCV[,'animal']/
  (constant30_herbs$VCV[,'animal']+constant30_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments30_shrubs <- constant_experiments30_final %>%
  filter(Growth_form == "Shrub")

constant30_shrubs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "multinomial2", pedigree = meta_constant30_tree, 
                                        prior = priors, data = constant_experiments30_shrubs,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant30_shrubs)
n_distinct(constant_experiments30_shrubs$animal)
nrow(constant_experiments30_shrubs)

lambda <- constant30_shrubs$VCV[,'animal']/
  (constant30_shrubs$VCV[,'animal']+constant30_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments30_restricted <- constant_experiments30_final %>%
  filter(Distribution == "Restricted")

constant30_restricted <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "multinomial2", pedigree = meta_constant30_tree, 
                                            prior = priors, data = constant_experiments30_restricted,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant30_restricted)
n_distinct(constant_experiments30_shrubs$animal)
nrow(constant_experiments30_restricted)

lambda <- constant30_restricted$VCV[,'animal']/
  (constant30_restricted$VCV[,'animal']+constant30_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments30_widespread <- constant_experiments30_final %>%
  filter(Distribution == "Widespread")

constant30_widespread <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "multinomial2", pedigree = meta_constant30_tree, 
                                            prior = priors, data = constant_experiments30_widespread,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant30_widespread)
n_distinct(constant_experiments30_widespread$animal)
nrow(constant_experiments30_widespread)

lambda <- constant30_widespread$VCV[,'animal']/
  (constant30_widespread$VCV[,'animal']+constant30_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments30_mesic <- constant_experiments30_final %>%
  filter(Microhabitat == "Mesic")

constant30_mesic <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant30_tree, 
                                       prior = priors, data = constant_experiments30_mesic,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant30_mesic)
n_distinct(constant_experiments30_mesic$animal)
nrow(constant_experiments30_mesic)

lambda <- constant30_mesic$VCV[,'animal']/
  (constant30_mesic$VCV[,'animal']+constant30_mesic$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments30_mesic_xeric <- constant_experiments30_final %>%
  filter(Microhabitat == "Mesic/Xeric")

constant30_mesic_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "multinomial2", pedigree = meta_constant30_tree, 
                                             prior = priors, data = constant_experiments30_mesic_xeric,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant30_mesic_xeric)
n_distinct(constant_experiments30_mesic_xeric$animal)
nrow(constant_experiments30_mesic_xeric)

lambda <- constant30_mesic_xeric$VCV[,'animal']/
  (constant30_mesic_xeric$VCV[,'animal']+constant30_mesic_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments30_xeric <- constant_experiments30_final %>%
  filter(Microhabitat == "Xeric")

constant30_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant30_tree, 
                                       prior = priors, data = constant_experiments30_xeric,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant30_xeric)
n_distinct(constant_experiments30_xeric$animal)
nrow(constant_experiments30_xeric)

lambda <- constant30_xeric$VCV[,'animal']/
  (constant30_xeric$VCV[,'animal']+constant30_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments30_mass <- constant_experiments30_final %>%
  filter(!is.na(Dry_mass))

constant30_mass <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition * scale(Dry_mass), 
                                      random = ~ animal + ID + n:ID, 
                                      family = "multinomial2", pedigree = meta_constant30_tree, 
                                      prior = priors, data = constant_experiments30_mass,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant30_mass)
n_distinct(constant_experiments30_mass$animal)
nrow(constant_experiments30_mass)

lambda <- constant30_mass$VCV[,'animal']/
  (constant30_mass$VCV[,'animal']+constant30_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# t50

constant_experiments_time30 <- constant_experiments_time %>%
  filter(Temperature == "30" | Temperature == "25")

exclude_constant_experiments_time30 <- constant_experiments_time %>%
  filter(Temperature == "30" | Temperature == "25") %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25" ~ "Control",
                                 Temperature == "30" ~ "Treatment")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 | Treatment < 0.1)) %>%
  select(1:3)

constant_experiments_time30_2 <- constant_experiments_time30 %>%
  anti_join(exclude_constant_experiments_time30) %>%
  mutate(D0 = 0,
         n = as.factor(row_number())) %>%
  select(n, ID:nSeeds, D0, D1:Species_acceptedLCVP)

constant_experiments_time30_t50 <- germination.indices(as.data.frame(constant_experiments_time30_2),
                                                       total.seeds.col = "nSeeds",
                                                       counts.intervals.cols = intervals30,
                                                       intervals = 1:31,
                                                       t50 = TRUE,
                                                       EmergenceRateIndex = FALSE,
                                                       PeakGermPercent = FALSE,
                                                       PeakGermTime = FALSE,
                                                       GermValue = FALSE) %>%
  dplyr::select(n:Temperature, Growth_form:Species_acceptedLCVP, t50_Farooq)

write.csv(constant_experiments_time30_t50, file = "constant_experiments_time30_missing.csv")

constant_experiments_time30_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/constant_experiments_time30_final.csv")) %>%
  filter(!is.na(Species_acceptedLCVP)) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Condition = case_when(Condition == "30" ~ 1,
                               Condition == "25" ~ 0),
         Condition = as.factor(Condition))

meta_constant_time30_names <- constant_experiments_time30_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_constant_time30_names <- names_for_pruning(meta_constant_time30_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_constant_time30_tree <- prune.sample(meta_constant_time30_names, phylo = full_tree)
meta_constant_time30_tree$node.label <- NULL

constant_time30_all <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                          random = ~ animal + ID + n:ID, 
                                          family = "gaussian", pedigree = meta_constant_time30_tree, 
                                          prior = priors, data = constant_experiments_time30_final,
                                          nitt = nite, thin = nthi, burnin = nbur, 
                                          verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                          pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time30_all)
n_distinct(constant_experiments_time30_final$animal)
nrow(constant_experiments_time30_final)

lambda <- constant_time30_all$VCV[,'animal']/
  (constant_time30_all$VCV[,'animal']+constant_time30_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time30_herbs <- constant_experiments_time30_final %>%
  filter(Growth_form == "Herb")

constant_time30_herbs <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time30_tree, 
                                            prior = priors, data = constant_experiments_time30_herbs,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time30_herbs)
n_distinct(constant_experiments_time30_herbs$animal)
nrow(constant_experiments_time30_herbs)

lambda <- constant_time20_herbs$VCV[,'animal']/
  (constant_time20_herbs$VCV[,'animal']+constant_time20_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time30_shrubs <- constant_experiments_time30_final %>%
  filter(Growth_form == "Shrub")

constant_time30_shrubs <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "gaussian", pedigree = meta_constant_time30_tree, 
                                             prior = priors, data = constant_experiments_time30_shrubs,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time30_shrubs)
n_distinct(constant_experiments_time30_shrubs$animal)
nrow(constant_experiments_time30_shrubs)

lambda <- constant_time30_shrubs$VCV[,'animal']/
  (constant_time30_shrubs$VCV[,'animal']+constant_time30_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time30_restricted <- constant_experiments_time30_final %>%
  filter(Distribution == "Restricted")

constant_time30_restricted <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                 random = ~ animal + ID + n:ID, 
                                                 family = "gaussian", pedigree = meta_constant_time30_tree, 
                                                 prior = priors, data = constant_experiments_time30_restricted,
                                                 nitt = nite, thin = nthi, burnin = nbur, 
                                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                 pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time30_restricted)
n_distinct(constant_experiments_time30_restricted$animal)
nrow(constant_experiments_time30_restricted)

lambda <- constant_time30_restricted$VCV[,'animal']/
  (constant_time30_restricted$VCV[,'animal']+constant_time30_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time30_widespread <- constant_experiments_time30_final %>%
  filter(Distribution == "Widespread")

constant_time30_widespread <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                 random = ~ animal + ID + n:ID, 
                                                 family = "gaussian", pedigree = meta_constant_time30_tree, 
                                                 prior = priors, data = constant_experiments_time30_widespread,
                                                 nitt = nite, thin = nthi, burnin = nbur, 
                                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                 pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time30_widespread)
n_distinct(constant_experiments_time30_widespread$animal)
nrow(constant_experiments_time30_widespread)

lambda <- constant_time30_widespread$VCV[,'animal']/
  (constant_time30_widespread$VCV[,'animal']+constant_time30_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time30_mesic <- constant_experiments_time30_final %>%
  filter(Microhabitat == "Mesic")

constant_time30_mesic <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time30_tree,
                                            prior = priors, data = constant_experiments_time30_mesic,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time30_mesic)
n_distinct(constant_experiments_time30_mesic$animal)
nrow(constant_experiments_time30_mesic)

lambda <- constant_time30_mesic$VCV[,'animal']/
  (constant_time30_mesic$VCV[,'animal']+constant_time30_mesic$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time30_mesic_xeric <- constant_experiments_time30_final %>%
  filter(Microhabitat == "Mesic/Xeric")

constant_time30_mesic_xeric <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                  random = ~ animal + ID + n:ID, 
                                                  family = "gaussian", pedigree = meta_constant_time30_tree,
                                                  prior = priors, data = constant_experiments_time30_mesic_xeric,
                                                  nitt = nite, thin = nthi, burnin = nbur, 
                                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time30_mesic_xeric)
n_distinct(constant_experiments_time30_mesic_xeric$animal)
nrow(constant_experiments_time30_mesic_xeric)

lambda <- constant_time30_mesic_xeric$VCV[,'animal']/
  (constant_time30_mesic_xeric$VCV[,'animal']+constant_time30_mesic_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time30_xeric <- constant_experiments_time30_final %>%
  filter(Microhabitat == "Xeric")

constant_time30_xeric <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time30_tree, 
                                            prior = priors, data = constant_experiments_time30_xeric,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time30_xeric)
n_distinct(constant_experiments_time30_xeric$animal)
nrow(constant_experiments_time30_xeric)

lambda <- constant_time30_xeric$VCV[,'animal']/
  (constant_time30_xeric$VCV[,'animal']+constant_time30_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time30_mass <- constant_experiments_time30_final %>%
  filter(!is.na(Dry_mass))

constant_time30_mass <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition * scale(Dry_mass), 
                                           random = ~ animal + ID + n:ID, 
                                           family = "gaussian", pedigree = meta_constant_time30_tree, 
                                           prior = priors, data = constant_experiments_time30_mass,
                                           nitt = nite, thin = nthi, burnin = nbur, 
                                           verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                           pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time30_mass)
n_distinct(constant_experiments_time30_mass$animal)
nrow(constant_experiments_time30_mass)

lambda <- constant_time30_mass$VCV[,'animal']/
  (constant_time30_mass$VCV[,'animal']+constant_time30_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

#### 7. 35 °C ####

constant_experiments35 <- constant_experiments %>%
  filter(Temperature == "35" | Temperature == "25")

exclude_constant_experiments35 <- constant_experiments %>%
  filter(Temperature == "35" | Temperature == "25") %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25" ~ "Control",
                                 Temperature == "35" ~ "Treatment")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 & Treatment < 0.1)) %>%
  select(1:3)

constant_experiments35_2 <- constant_experiments35 %>%
  anti_join(exclude_constant_experiments35)

# Creating the dataset for adding missing seed mass values

write.csv(constant_experiments35_2, file = "constant_experiments35_missing.csv")

# Loading the dataset with the inputed mass data.

constant_experiments35_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/constant_experiments35_final.csv")) %>%
  filter(!is.na(Species_acceptedLCVP)) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Condition = case_when(Condition == "35" ~ 1,
                               Condition == "25" ~ 0),
         Condition = as.factor(Condition))

meta_constant35_names <- constant_experiments35_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_constant35_names <- names_for_pruning(meta_constant35_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_constant35_tree <- prune.sample(meta_constant35_names, phylo = full_tree)
meta_constant35_tree$node.label <- NULL

constant35_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                     random = ~ animal + ID + n:ID, 
                                     family = "multinomial2", pedigree = meta_constant35_tree, 
                                     prior = priors, data = constant_experiments35_final,
                                     nitt = nite, thin = nthi, burnin = nbur, 
                                     verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                     pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant35_all)
n_distinct(constant_experiments35_final$animal)
nrow(constant_experiments35_final)

lambda <- constant35_all$VCV[,'animal']/
  (constant35_all$VCV[,'animal']+constant35_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments35_herbs <- constant_experiments35_final %>%
  filter(Growth_form == "Herb")

constant35_herbs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant35_tree, 
                                       prior = priors, data = constant_experiments35_herbs,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant35_herbs)
n_distinct(constant_experiments35_herbs$animal)
nrow(constant_experiments35_herbs)

lambda <- constant35_herbs$VCV[,'animal']/
  (constant35_herbs$VCV[,'animal']+constant35_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments35_shrubs <- constant_experiments35_final %>%
  filter(Growth_form == "Shrub")

constant35_shrubs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "multinomial2", pedigree = meta_constant35_tree, 
                                        prior = priors, data = constant_experiments35_shrubs,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant35_shrubs)
n_distinct(constant_experiments35_shrubs$animal)
nrow(constant_experiments35_shrubs)

lambda <- constant35_shrubs$VCV[,'animal']/
  (constant35_shrubs$VCV[,'animal']+constant35_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments35_restricted <- constant_experiments35_final %>%
  filter(Distribution == "Restricted")

constant35_restricted <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "multinomial2", pedigree = meta_constant35_tree, 
                                            prior = priors, data = constant_experiments35_restricted,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant35_restricted)
n_distinct(constant_experiments35_restricted$animal)
nrow(constant_experiments35_restricted)

lambda <- constant35_restricted$VCV[,'animal']/
  (constant35_restricted$VCV[,'animal']+constant35_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments35_widespread <- constant_experiments35_final %>%
  filter(Distribution == "Widespread")

constant35_widespread <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "multinomial2", pedigree = meta_constant35_tree, 
                                            prior = priors, data = constant_experiments35_widespread,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant35_widespread)
n_distinct(constant_experiments35_widespread$animal)
nrow(constant_experiments35_widespread)

lambda <- constant35_widespread$VCV[,'animal']/
  (constant35_widespread$VCV[,'animal']+constant35_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments35_mesic <- constant_experiments35_final %>%
  filter(Microhabitat == "Mesic")

constant35_mesic <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant35_tree, 
                                       prior = priors, data = constant_experiments35_mesic,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant35_mesic)
n_distinct(constant_experiments35_mesic$animal)
nrow(constant_experiments35_mesic)

lambda <- constant35_mesic$VCV[,'animal']/
  (constant35_mesic$VCV[,'animal']+constant35_mesic$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments35_mesic_xeric <- constant_experiments35_final %>%
  filter(Microhabitat == "Mesic/Xeric")

constant35_mesic_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "multinomial2", pedigree = meta_constant35_tree, 
                                             prior = priors, data = constant_experiments35_mesic_xeric,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant35_mesic_xeric)
n_distinct(constant_experiments35_mesic_xeric$animal)
nrow(constant_experiments35_mesic_xeric)

lambda <- constant35_mesic_xeric$VCV[,'animal']/
  (constant35_mesic_xeric$VCV[,'animal']+constant35_mesic_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments35_xeric <- constant_experiments35_final %>%
  filter(Microhabitat == "Xeric")

constant35_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = meta_constant35_tree, 
                                       prior = priors, data = constant_experiments35_xeric,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant35_xeric)
n_distinct(constant_experiments35_xeric$animal)
nrow(constant_experiments35_xeric)

lambda <- constant35_xeric$VCV[,'animal']/
  (constant35_xeric$VCV[,'animal']+constant35_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments35_mass <- constant_experiments35_final %>%
  filter(!is.na(Dry_mass))

constant35_mass <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition * scale(Dry_mass), 
                                      random = ~ animal + ID + n:ID, 
                                      family = "multinomial2", pedigree = meta_constant35_tree, 
                                      prior = priors, data = constant_experiments35_mass,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant35_mass)
n_distinct(constant_experiments35_mass$animal)
nrow(constant_experiments35_mass)

lambda <- constant35_mass$VCV[,'animal']/
  (constant35_mass$VCV[,'animal']+constant35_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# t50

constant_experiments_time35 <- constant_experiments_time %>%
  filter(Temperature == "35" | Temperature == "25")

exclude_constant_experiments_time35 <- constant_experiments_time %>%
  filter(Temperature == "35" | Temperature == "25") %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25" ~ "Control",
                                 Temperature == "35" ~ "Treatment")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 | Treatment < 0.1)) %>%
  select(1:3)

constant_experiments_time35_2 <- constant_experiments_time35 %>%
  anti_join(exclude_constant_experiments_time35) %>%
  mutate(D0 = 0,
         n = as.factor(row_number())) %>%
  select(n, ID:nSeeds, D0, D1:Species_acceptedLCVP)

constant_experiments_time35_t50 <- germination.indices(as.data.frame(constant_experiments_time35_2),
                                                       total.seeds.col = "nSeeds",
                                                       counts.intervals.cols = intervals30,
                                                       intervals = 1:31,
                                                       t50 = TRUE,
                                                       EmergenceRateIndex = FALSE,
                                                       PeakGermPercent = FALSE,
                                                       PeakGermTime = FALSE,
                                                       GermValue = FALSE) %>%
  dplyr::select(n:Temperature, Growth_form:Species_acceptedLCVP, t50_Farooq)

write.csv(constant_experiments_time35_t50, file = "constant_experiments_time35_missing.csv")

constant_experiments_time35_final <- read.csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/constant_experiments_time35_final.csv")) %>%
  filter(!is.na(Species_acceptedLCVP)) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Condition = case_when(Condition == "35" ~ 1,
                               Condition == "25" ~ 0),
         Condition = as.factor(Condition))

meta_constant_time35_names <- constant_experiments_time35_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_constant_time35_names <- names_for_pruning(meta_constant_time35_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_constant_time35_tree <- prune.sample(meta_constant_time35_names, phylo = full_tree)
meta_constant_time35_tree$node.label <- NULL

constant_time35_all <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                          random = ~ animal + ID + n:ID, 
                                          family = "gaussian", pedigree = meta_constant_time35_tree, 
                                          prior = priors, data = constant_experiments_time35_final,
                                          nitt = nite, thin = nthi, burnin = nbur, 
                                          verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                          pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time35_all)
n_distinct(constant_experiments_time35_final$animal)
nrow(constant_experiments_time35_final)

lambda <- constant_time35_all$VCV[,'animal']/
  (constant_time35_all$VCV[,'animal']+constant_time35_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time35_herbs <- constant_experiments_time35_final %>%
  filter(Growth_form == "Herb")

constant_time35_herbs <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time35_tree, 
                                            prior = priors, data = constant_experiments_time35_herbs,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time35_herbs)
n_distinct(constant_experiments_time35_herbs$animal)
nrow(constant_experiments_time35_herbs)

lambda <- constant_time20_herbs$VCV[,'animal']/
  (constant_time20_herbs$VCV[,'animal']+constant_time20_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time35_shrubs <- constant_experiments_time35_final %>%
  filter(Growth_form == "Shrub")

constant_time35_shrubs <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "gaussian", pedigree = meta_constant_time35_tree, 
                                             prior = priors, data = constant_experiments_time35_shrubs,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time35_shrubs)
n_distinct(constant_experiments_time35_shrubs$animal)
nrow(constant_experiments_time35_shrubs)

lambda <- constant_time35_shrubs$VCV[,'animal']/
  (constant_time35_shrubs$VCV[,'animal']+constant_time35_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time35_restricted <- constant_experiments_time35_final %>%
  filter(Distribution == "Restricted")

constant_time35_restricted <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                 random = ~ animal + ID + n:ID, 
                                                 family = "gaussian", pedigree = meta_constant_time35_tree, 
                                                 prior = priors, data = constant_experiments_time35_restricted,
                                                 nitt = nite, thin = nthi, burnin = nbur, 
                                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                 pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time35_restricted)
n_distinct(constant_experiments_time35_restricted$animal)
nrow(constant_experiments_time35_restricted)

lambda <- constant_time35_restricted$VCV[,'animal']/
  (constant_time35_restricted$VCV[,'animal']+constant_time35_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time35_widespread <- constant_experiments_time35_final %>%
  filter(Distribution == "Widespread")

constant_time35_widespread <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                 random = ~ animal + ID + n:ID, 
                                                 family = "gaussian", pedigree = meta_constant_time35_tree, 
                                                 prior = priors, data = constant_experiments_time35_widespread,
                                                 nitt = nite, thin = nthi, burnin = nbur, 
                                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                 pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time35_widespread)
n_distinct(constant_experiments_time35_widespread$animal)
nrow(constant_experiments_time35_widespread)

lambda <- constant_time35_widespread$VCV[,'animal']/
  (constant_time35_widespread$VCV[,'animal']+constant_time35_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time35_mesic <- constant_experiments_time35_final %>%
  filter(Microhabitat == "Mesic")

constant_time35_mesic <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time35_tree,
                                            prior = priors, data = constant_experiments_time35_mesic,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time35_mesic)
n_distinct(constant_experiments_time35_mesic$animal)
nrow(constant_experiments_time35_mesic)

lambda <- constant_time35_mesic$VCV[,'animal']/
  (constant_time35_mesic$VCV[,'animal']+constant_time35_mesic$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time35_mesic_xeric <- constant_experiments_time35_final %>%
  filter(Microhabitat == "Mesic/Xeric")

constant_time35_mesic_xeric <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                                  random = ~ animal + ID + n:ID, 
                                                  family = "gaussian", pedigree = meta_constant_time35_tree,
                                                  prior = priors, data = constant_experiments_time35_mesic_xeric,
                                                  nitt = nite, thin = nthi, burnin = nbur, 
                                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time35_mesic_xeric)
n_distinct(constant_experiments_time35_mesic_xeric$animal)
nrow(constant_experiments_time35_mesic_xeric)

lambda <- constant_time35_mesic_xeric$VCV[,'animal']/
  (constant_time35_mesic_xeric$VCV[,'animal']+constant_time35_mesic_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time35_xeric <- constant_experiments_time35_final %>%
  filter(Microhabitat == "Xeric")

constant_time35_xeric <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = meta_constant_time35_tree, 
                                            prior = priors, data = constant_experiments_time35_xeric,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time35_xeric)
n_distinct(constant_experiments35_xeric$animal)
nrow(constant_experiments_time35_xeric)

lambda <- constant_time35_xeric$VCV[,'animal']/
  (constant_time35_xeric$VCV[,'animal']+constant_time35_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments_time35_mass <- constant_experiments_time35_final %>%
  filter(!is.na(Dry_mass))

constant_time35_mass <- MCMCglmm::MCMCglmm(t50_Farooq ~ Condition * scale(Dry_mass), 
                                           random = ~ animal + ID + n:ID, 
                                           family = "gaussian", pedigree = meta_constant_time35_tree, 
                                           prior = priors, data = constant_experiments_time35_mass,
                                           nitt = nite, thin = nthi, burnin = nbur, 
                                           verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                           pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant_time35_mass)
n_distinct(constant_experiments_time35_mass$animal)
nrow(constant_experiments_time35_mass)

lambda <- constant_time35_mass$VCV[,'animal']/
  (constant_time35_mass$VCV[,'animal']+constant_time35_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

#### 8. 40 °C ####

constant_experiments40 <- constant_experiments %>%
  filter(Temperature == "40" | Temperature == "25")

exclude_constant_experiments40 <- constant_experiments %>%
  filter(Temperature == "40" | Temperature == "25") %>%
  dplyr::select(1:6) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Temperature = case_when(Temperature == "25" ~ "Control",
                                 Temperature == "40" ~ "Treatment")) %>%
  group_by(ID, Species_reported, Species_accepted, Temperature) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Temperature, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 & Treatment < 0.1)) %>%
  select(1:3)

constant_experiments40_2 <- constant_experiments40 %>%
  anti_join(exclude_constant_experiments40)

# Creating the dataset for adding missing seed mass values

write.csv(constant_experiments40_2, file = "constant_experiments40_missing.csv")

# Loading the dataset with the inputed mass data.

constant_experiments40_final <- read_csv(str_c(git_find(),"/Analyses/Meta_Analyses/Temperature/constant_experiments40_final.csv")) %>%
  filter(!is.na(Species_acceptedLCVP)) %>%
  rename(Condition = Temperature,
         animal = Species_acceptedLCVP) %>%
  mutate(n = as.factor(row_number()),
         Condition = case_when(Condition == "40" ~ 1,
                               Condition == "25" ~ 0),
         Condition = as.factor(Condition))

constant_experiments40_final <- as.data.frame(constant_experiments40_final)

meta_constant40_names <- constant_experiments40_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

meta_constant40_names <-  names_for_pruning(meta_constant40_names)

# Now, we prune the phylogeny in order to keep only these names.

meta_constant40_tree <- prune.sample(meta_constant40_names, phylo = full_tree)
meta_constant40_tree$node.label <- NULL

constant40_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                     random = ~ animal + ID + n:ID, 
                                     family = "multinomial2", pedigree = meta_constant40_tree, 
                                     prior = priors, data = constant_experiments40_final,
                                     nitt = nite, thin = nthi, burnin = nbur, 
                                     verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                     pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant40_all)
n_distinct(constant_experiments40_final$animal)
nrow(constant_experiments40_final)

lambda <- constant40_all$VCV[,'animal']/
  (constant40_all$VCV[,'animal']+constant40_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments40_herbs <- constant_experiments40_final %>%
  filter(Growth_form == "Herb")

constant40_herbs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + n, 
                                       family = "multinomial2", pedigree = meta_constant40_tree, 
                                       prior = priors2, data = constant_experiments40_herbs,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant40_herbs)
n_distinct(constant_experiments40_herbs$animal)
nrow(constant_experiments40_final)

lambda <- constant40_herbs$VCV[,'animal']/
  (constant40_herbs$VCV[,'animal']+constant40_herbs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments40_shrubs <- constant_experiments40_final %>%
  filter(Growth_form == "Shrub")

constant40_shrubs <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "multinomial2", pedigree = meta_constant40_tree, 
                                        prior = priors, data = constant_experiments40_shrubs,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant40_shrubs)
n_distinct(constant_experiments40_shrubs$animal)
nrow(constant_experiments40_shrubs)

lambda <- constant40_shrubs$VCV[,'animal']/
  (constant40_shrubs$VCV[,'animal']+constant40_shrubs$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments40_restricted <- constant_experiments40_final %>%
  filter(Distribution == "Restricted")

constant40_restricted <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + n, 
                                            family = "multinomial2", pedigree = meta_constant40_tree, 
                                            prior = priors2, data = constant_experiments40_restricted,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant40_restricted)
n_distinct(constant_experiments40_shrubs$animal)
nrow(constant_experiments40_restricted)

lambda <- constant40_restricted$VCV[,'animal']/
  (constant40_restricted$VCV[,'animal']+constant40_restricted$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments40_widespread <- constant_experiments40_final %>%
  filter(Distribution == "Widespread")

constant40_widespread <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "multinomial2", pedigree = meta_constant40_tree, 
                                            prior = priors, data = constant_experiments40_widespread,
                                            nitt = nite, thin = nthi, burnin = nbur, 
                                            verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                            pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant40_widespread)
n_distinct(constant_experiments40_widespread$animal)
nrow(constant_experiments40_widespread)

lambda <- constant40_widespread$VCV[,'animal']/
  (constant40_widespread$VCV[,'animal']+constant40_widespread$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments40_mesic <- constant_experiments40_final %>%
  filter(Microhabitat == "Mesic")

constant40_mesic <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + n, 
                                       family = "multinomial2", pedigree = meta_constant40_tree, 
                                       prior = priors2, data = constant_experiments40_mesic,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant40_mesic)
n_distinct(constant_experiments40_mesic$animal)
nrow(constant_experiments40_mesic)

lambda <- constant40_mesic$VCV[,'animal']/
  (constant40_mesic$VCV[,'animal']+constant40_mesic$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments40_mesic_xeric <- constant_experiments40_final %>%
  filter(Microhabitat == "Mesic/Xeric")

constant40_mesic_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "multinomial2", pedigree = meta_constant40_tree, 
                                             prior = priors, data = constant_experiments40_mesic_xeric,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant40_mesic_xeric)
n_distinct(constant_experiments40_mesic_xeric$animal)
nrow(constant_experiments40_mesic_xeric)

lambda <- constant40_mesic_xeric$VCV[,'animal']/
  (constant40_mesic_xeric$VCV[,'animal']+constant40_mesic_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments40_xeric <- constant_experiments40_final %>%
  filter(Microhabitat == "Xeric")

constant40_xeric <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition, 
                                       random = ~ animal + n, 
                                       family = "multinomial2", pedigree = meta_constant40_tree, 
                                       prior = priors2, data = constant_experiments40_xeric,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant40_xeric)
n_distinct(constant_experiments40_xeric$animal)
nrow(constant_experiments40_xeric)

lambda <- constant40_xeric$VCV[,'animal']/
  (constant40_xeric$VCV[,'animal']+constant40_xeric$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

constant_experiments40_mass <- constant_experiments40_final %>%
  filter(!is.na(Dry_mass))

constant40_mass <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Condition * scale(Dry_mass), 
                                      random = ~ animal + ID + n:ID, 
                                      family = "multinomial2", pedigree = meta_constant40_tree, 
                                      prior = priors, data = constant_experiments40_mass,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(constant40_mass)
n_distinct(constant_experiments40_mass$animal)
nrow(constant_experiments40_mass)

lambda <- constant40_mass$VCV[,'animal']/
  (constant40_mass$VCV[,'animal']+constant40_mass$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

