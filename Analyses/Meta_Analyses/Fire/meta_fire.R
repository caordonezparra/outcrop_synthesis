# The following code was designed to assess the effect of fire-related cues on germination percentage and median germination time by performing meta-analyses through fitting binomial or Gaussian phylogenetic generalized mixed models with Bayesian estimation using Markov chain Monte Carlo (MCMC) method as implemented in the MCMCglmm R-package.

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: April 28th, 2024.

#### 0. Loading packages (assuming they are already installed) ####

library(readr) # v. 2.1.5. for database loading.
library(dplyr) # v. 1.1.4 for database manipulation.
library(tidyr) # v. 1.3.0. (also) for database manipulation
library(tibble) # v. 3.2.1 (also) for database manipulation
library(forcats) # v. 1.0.0 for working with categorical variables.
library(purrr) # v. 1.0.2 for reiteratively applicating functions.
library(picante) # v. 1.8.2 for loading phylogenetic trees
library(germinationmetrics) # v. 0.1.8.9000 to calculate t50
library(MCMCglmm) # v. 2.35 for generalized mixed models with Bayesian estimation using MCMC method.

#### 1. Setting working directory

setwd(str_c(git_find(),"/Analyses/Meta_Analyses/Fire"))

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

#### 2. Getting the data from studies that assessed the effect of fire ####

# The following code looks to create a data set of the germination experiments that tested the effect of fire-related cues. The control treatment in all cases were untreated seeds (i.e., seeds not submitted to heat shock nor smoke).

# Heat

# First we're getting the list of studies assessing the effects of heat shocks. Therefore, we are filtering the "topics" database looking for  studies where "Minor" is "Heat".

heat <- topics %>%
  filter(Minor == "Heat")

# Then, we are getting the IDs of each of these studies.

heat_studies <- unique(heat$ID)
heat_studies

# Now, we are filtering the germination database, looking for:

heat_experiments <- germination %>%
  filter(ID %in% heat_studies, # Studies with whose IDs are in "heat_studies",
         Photoperiod != 0,
         !is.na(HeatShock_Temperature), # that they evaluated the effect of heat
         is.na(Smoke) | Smoke == "Control") %>% # shocks alone (i.e., not in combination with smoke)
  select(1:3, HeatShock_Temperature, HeatShock_Time, Smoke, GermSeeds, nSeeds) %>%
  left_join(traits) %>%
  select(1:8, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits info.
  left_join(lcvp_names_checked) %>%
  filter(HeatShock_Temperature == "Control" | HeatShock_Temperature %in% c(80, 100, 200),
         Growth_form == "Herb" | Growth_form == "Shrub",
         Dormancy != "NC") %>%
  unite("HeatShock_Treatment", HeatShock_Temperature, HeatShock_Time, sep = "*", remove = FALSE) %>%
  mutate(HeatShock_Treatment = case_when(HeatShock_Treatment == "Control*Control" ~ "Control",
                                         HeatShock_Treatment != "Control*Control" ~ "Treatment"))

exclude_heat_experiments <- heat_experiments %>%
  select(1:6, HeatShock_Treatment, GermSeeds, nSeeds) %>%
  mutate(GermProp = GermSeeds/nSeeds) %>%
  group_by(ID, Species_reported, Species_accepted, HeatShock_Treatment) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = HeatShock_Treatment, values_from = meanGermProp)%>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 & Treatment < 0.1)) %>%
  select(1:3)

heat_experiments_2 <- heat_experiments %>%
  anti_join(exclude_heat_experiments) %>%
  unite("HeatShock_Treatment", HeatShock_Temperature, HeatShock_Time, sep = "*", remove = FALSE) %>%
  mutate(Treatment_type = "Single",
         Treatment = "Heat",
         HeatShock_Treatment = case_when(HeatShock_Treatment == "Control*Control" ~ "Control",
                                                HeatShock_Treatment != "Control*Control" ~ HeatShock_Treatment),
         Treatment_application = case_when(HeatShock_Treatment == "Control" ~ 0,
                                           .default = 1)) %>%
  select(1:3, Treatment_type, Treatment, Treatment_application, HeatShock_Treatment, Smoke:Species_acceptedLCVP)

### Smoke

smoke <- topics %>%
  filter(Minor == "Smoke")

# Then, we are getting the IDs of each of these studies.

smoke_studies <- unique(smoke$ID)

# Now, we are filtering the germination database: 

smoke_experiments <- germination %>%
  filter(ID %in% smoke_studies, # Studies with IDs in "smoke_studies",
         Photoperiod != 0,
         Smoke == "Smoke water" | Smoke == "Control",
         is.na(HeatShock_Temperature) | HeatShock_Temperature == "Control") %>%  # Assessed only the effect of smoke
  select(1:3, HeatShock_Temperature, HeatShock_Time, Smoke, GermSeeds, nSeeds) %>%
  left_join(traits) %>%
  select(1:8, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits info.
  left_join(lcvp_names_checked) %>%
  filter(Growth_form == "Herb" | Growth_form == "Shrub",
         Dormancy != "NC") %>%
  mutate(HeatShock_Treatment = NA) %>%
  select(1:3, HeatShock_Treatment, 4:14)
  
exclude_smoke_experiments <- smoke_experiments %>%
  select(1:3, Smoke, GermSeeds, nSeeds) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Smoke = case_when(Smoke == "Smoke water" ~ "Treatment",
                           .default = "Control")) %>%
  group_by(ID, Species_reported, Species_accepted, Smoke) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Smoke, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 & Treatment < 0.1)) %>%
  select(1:3)

smoke_experiments_2 <- smoke_experiments %>%
  anti_join(exclude_smoke_experiments) %>%
  mutate(Treatment_type = "Single",
         Treatment = "Smoke",
         Treatment_application = case_when(Smoke == "Control" ~ 0,
                                           .default = 1)) %>%
  select(1:3, Treatment_type, Treatment, Treatment_application, HeatShock_Treatment, Smoke:Species_acceptedLCVP)

# Heat shock and smoke combined

heat_smoke_studies <- intersect(heat_studies, smoke_studies)

heat_smoke_experiments <- germination %>%
  filter(ID %in% heat_smoke_studies,
         Photoperiod != 0,
         Smoke == "Smoke water" | Smoke == "Control",
         !is.na(HeatShock_Temperature)) %>%
  select(1:3, HeatShock_Temperature, HeatShock_Time, Smoke, GermSeeds, nSeeds) %>%
  left_join(traits) %>%
  select(1:8, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits info.
  left_join(lcvp_names_checked) %>%
  filter(HeatShock_Temperature == "Control" | HeatShock_Temperature %in% c(80, 100, 200),
         Growth_form == "Herb" | Growth_form == "Shrub",
         Dormancy != "NC") %>%
  unite("Treatment", HeatShock_Temperature, HeatShock_Time, Smoke, sep = "*", remove = FALSE) %>%
  mutate(Treatment = case_when(Treatment == "Control*Control*Control" ~ "Control",
                               Treatment != "Control*Control*Control" ~ "Treatment"))

exclude_heat_smoke_experiments <- heat_smoke_experiments %>%
  select(1:3, Treatment, GermSeeds, nSeeds) %>%
  mutate(GermProp = GermSeeds/nSeeds) %>%
  group_by(ID, Species_reported, Species_accepted, Treatment) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Treatment, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 & Treatment < 0.1)) %>%
  select(1:3)

heat_smoke_experiments_2 <- heat_smoke_experiments %>%
  anti_join(exclude_heat_smoke_experiments) %>%
  unite("HeatShock_Treatment", HeatShock_Temperature, HeatShock_Time, sep = "*", remove = FALSE) %>%
  mutate(Treatment_type = "Combined",
         Treatment = "Combined",
         Treatment_application = case_when(HeatShock_Treatment == "Control*Control" ~ 0,
                                           .default = 1)) %>%
  select(1:3, Treatment_type, Treatment, Treatment_application, HeatShock_Treatment, Smoke:Species_acceptedLCVP)

fire_cues_list <- list(heat_experiments_2, smoke_experiments_2, heat_smoke_experiments_2)

fire_cues <- fire_cues_list %>%
  reduce(full_join) %>%
  filter(!is.na(Species_acceptedLCVP))
  
# We are going to export this dataset to include the missing data of seed mass

write.csv(fire_cues, "fire_cues_missing.csv")

fire_cues_final <- as.data.frame(read_csv(str_c(git_find(),"/Analyses/Meta_Analyses/Fire/fire_cues_final.csv"))) %>%
  mutate(ID = as.factor(ID),
         Treatment_application = as.factor(Treatment_application),
         n = row_number(),
         n = as.factor(n)) %>%
  rename(animal = Species_acceptedLCVP)

n_distinct(fire_cues_final$animal)
nrow(fire_cues_final)

fire_cues_names <- fire_cues_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

fire_cues_names <- names_for_pruning(fire_cues_names)

# Now, we prune the phylogeny in order to keep only these names.

fire_cues_tree <- prune.sample(fire_cues_names, phylo = full_tree)
fire_cues_tree$node.label <- NULL

# Now, we will set the number of desired iterations based on De Villemereuil & Nakagawa (2014)

nite = 500000
nthi = 50
nbur = 50000

# Set priors based on https://doi.org/10.1017/S0960258517000332

priors <- list(R = list(V = 1, nu = 50), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))

#### The effect of fire-related cues all together ###

# We are going to test the overall effect of fire-related cues on germination.


fire_cues_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                     random = ~ animal + ID + n:ID, 
                                     family = "multinomial2", pedigree = fire_cues_tree, 
                                     prior = priors, data = fire_cues_final,
                                     nitt = nite, thin = nthi, burnin = nbur, 
                                     verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                     pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_all)
n_distinct(fire_cues_final$animal)

lambda <- fire_cues_all$VCV[,'animal']/
  (fire_cues_all$VCV[,'animal']+fire_cues_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_combined <- fire_cues_final %>%
  filter(Treatment_type == "Combined")

n_distinct(fire_cues_combined$animal)
nrow(fire_cues_combined)

fire_cues_combinedt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                    random = ~ animal + n, 
                                    family = "multinomial2", pedigree = fire_cues_tree, 
                                    prior = priors2, data = fire_cues_combined,
                                    nitt = nite, thin = nthi, burnin = nbur, 
                                    verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                    pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_combinedt)
n_distinct(fire_cues_combined$animal)

lambda <- fire_cues_combinedt$VCV[,'animal']/
  (fire_cues_combinedt$VCV[,'animal']+fire_cues_combinedt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_single <- fire_cues_final %>%
  filter(Treatment_type == "Single")

n_distinct(fire_cues_single$animal)
nrow(fire_cues_single)

fire_cues_singlet <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                          random = ~ animal + ID + n:ID, 
                                          family = "multinomial2", pedigree = fire_cues_tree, 
                                          prior = priors, data = fire_cues_single,
                                          nitt = nite, thin = nthi, burnin = nbur, 
                                          verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                          pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_singlet)
n_distinct(fire_cues_single$animal)

lambda <- fire_cues_singlet$VCV[,'animal']/
  (fire_cues_singlet$VCV[,'animal']+fire_cues_singlet$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_heat <- fire_cues_final %>%
  filter(Treatment == "Heat") %>%
  mutate(HeatShock_Treatment = as.factor(HeatShock_Treatment))

n_distinct(fire_cues_heat$animal)
nrow(fire_cues_heat)

fire_cues_heatt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "multinomial2", pedigree = fire_cues_tree, 
                                        prior = priors, data = fire_cues_heat,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_heatt)
n_distinct(fire_cues_final$animal)

lambda <- fire_cues_heatt$VCV[,'animal']/
  (fire_cues_heatt$VCV[,'animal']+fire_cues_heatt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_smoke <- fire_cues_final %>%
  filter(Treatment == "Smoke")

n_distinct(fire_cues_smoke$animal)
nrow(fire_cues_smoke)

fire_cues_smoket <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                      random = ~ animal + ID + n:ID, 
                                      family = "multinomial2", pedigree = fire_cues_tree, 
                                      prior = priors, data = fire_cues_smoke,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_smoket)
n_distinct(fire_cues_smoke$animal)

lambda <- fire_cues_smoket$VCV[,'animal']/
  (fire_cues_smoket$VCV[,'animal']+fire_cues_smoket$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_dormant <- fire_cues_final %>%
  filter(Dormancy == "D")

n_distinct(fire_cues_dormant$animal)
nrow(fire_cues_dormant)

fire_cues_dormantt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = fire_cues_tree, 
                                       prior = priors, data = fire_cues_dormant,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_dormantt)
n_distinct(fire_cues_dormant$animal)

lambda <- fire_cues_dormantt$VCV[,'animal']/
  (fire_cues_dormantt$VCV[,'animal']+fire_cues_dormantt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_nondormant <- fire_cues_final %>%
  filter(Dormancy == "ND")

n_distinct(fire_cues_nondormant$animal)
nrow(fire_cues_nondormant)

fire_cues_non_dormantt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                         random = ~ animal + ID + n:ID, 
                                         family = "multinomial2", pedigree = fire_cues_tree, 
                                         prior = priors, data = fire_cues_nondormant,
                                         nitt = nite, thin = nthi, burnin = nbur, 
                                         verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                         pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_non_dormantt)
n_distinct(fire_cues_nondormant$animal)

lambda <- fire_cues_non_dormantt$VCV[,'animal']/
  (fire_cues_non_dormantt$VCV[,'animal']+fire_cues_non_dormantt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_mass <- fire_cues_final %>%
  filter(!is.na(Dry_mass))

n_distinct(fire_cues_mass$animal)
nrow(fire_cues_mass)

fire_cues_masst <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application * scale(Dry_mass), 
                                             random = ~ animal + ID + n:ID, 
                                             family = "multinomial2", pedigree = fire_cues_tree, 
                                             prior = priors, data = fire_cues_mass,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_masst)
n_distinct(fire_cues_mass$animal)

lambda <- fire_cues_masst$VCV[,'animal']/
  (fire_cues_masst$VCV[,'animal']+fire_cues_masst$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Now we will test the effect of different heat shock temperatures

fire_cues_heat2 <- fire_cues_heat %>%
  filter(HeatShock_Treatment != "100*3") # We are excluding this treatment as it was only tested in two species in a      single study.

figure_cues_80_2 <- fire_cues_heat2 %>%
  filter(HeatShock_Treatment == "Control" | HeatShock_Treatment == "80*2")

n_distinct(figure_cues_80_2 $animal)
nrow(figure_cues_80_2)

fire_cues_80_2t <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "multinomial2", pedigree = fire_cues_tree, 
                                        prior = priors, data = figure_cues_80_2,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_80_2t)
n_distinct(figure_cues_80_2$animal)

lambda <- fire_cues_80_2t$VCV[,'animal']/
  (fire_cues_80_2t$VCV[,'animal']+fire_cues_80_2t$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

figure_cues_80_5 <- fire_cues_heat2 %>%
  filter(HeatShock_Treatment == "Control" | HeatShock_Treatment == "80*5")

n_distinct(figure_cues_80_5$animal)
nrow(figure_cues_80_5)

fire_cues_80_5t <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                      random = ~ animal + ID + n:ID, 
                                      family = "multinomial2", pedigree = fire_cues_tree, 
                                      prior = priors, data = figure_cues_80_5,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_80_5t)
n_distinct(figure_cues_80_5$animal)

lambda <- fire_cues_80_5t$VCV[,'animal']/
  (fire_cues_80_5t$VCV[,'animal']+fire_cues_80_5t$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

figure_cues_100_1 <- fire_cues_heat2 %>%
  filter(HeatShock_Treatment == "Control" | HeatShock_Treatment == "100*1")

n_distinct(figure_cues_100_1$animal)
nrow(figure_cues_100_1)

fire_cues_100_1t <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                      random = ~ animal + ID + n:ID, 
                                      family = "multinomial2", pedigree = fire_cues_tree, 
                                      prior = priors, data = figure_cues_100_1,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_100_1t)
n_distinct(figure_cues_100_1$animal)

lambda <- fire_cues_100_1t$VCV[,'animal']/
  (fire_cues_100_1t$VCV[,'animal']+fire_cues_100_1t$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

figure_cues_100_5 <- fire_cues_heat2 %>%
  filter(HeatShock_Treatment == "Control" | HeatShock_Treatment == "100*5")

n_distinct(figure_cues_100_5$animal)
nrow(figure_cues_100_5)

fire_cues_100_5t <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = fire_cues_tree, 
                                       prior = priors, data = figure_cues_100_5,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_100_5t)
n_distinct(figure_cues_100_5$animal)

lambda <- fire_cues_100_5t$VCV[,'animal']/
  (fire_cues_100_5t$VCV[,'animal']+fire_cues_100_5t$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

figure_cues_200_1 <- fire_cues_heat2 %>%
  filter(HeatShock_Treatment == "Control" | HeatShock_Treatment == "200*1")

n_distinct(figure_cues_200_1$animal)
nrow(figure_cues_200_1)

fire_cues_200_1t <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = fire_cues_tree, 
                                       prior = priors, data = figure_cues_200_1,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_200_1t)
n_distinct(figure_cues_200_1$animal)

lambda <- fire_cues_200_1t$VCV[,'animal']/
  (fire_cues_200_1t$VCV[,'animal']+fire_cues_200_1t$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Now, we'll be testing the effect of heat shocks on the different plant groups

fire_cues_heat3 <- fire_cues_heat %>%
  filter(HeatShock_Treatment != "200*1")

n_distinct(fire_cues_heat3$animal)
nrow(fire_cues_heat3)

heat_all <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                  random = ~ animal + ID + n:ID, 
                                  family = "multinomial2", pedigree = fire_cues_tree, 
                                  prior = priors, data = fire_cues_heat3,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_all)
n_distinct(fire_cues_heat3$animal)

lambda <- heat_all$VCV[,'animal']/
  (heat_all$VCV[,'animal']+heat_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_herbs <- fire_cues_heat3 %>%
  filter(Growth_form == "Herb")

n_distinct(heat_herbs$animal)
nrow(heat_herbs)

heat_herbst <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = fire_cues_tree, 
                                       prior = priors, data = heat_herbs,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_herbst)
n_distinct(heat_herbs$animal)

lambda <- heat_herbst$VCV[,'animal']/
  (heat_herbst$VCV[,'animal']+heat_herbst$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_shrub <- fire_cues_heat3 %>%
  filter(Growth_form == "Shrub")

n_distinct(heat_shrub$animal)
nrow(heat_shrub)

heat_shrubt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                  random = ~ animal + ID + n:ID, 
                                  family = "multinomial2", pedigree = fire_cues_tree, 
                                  prior = priors, data = heat_shrub,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_shrubt)
n_distinct(heat_shrub$animal)

lambda <- heat_shrubt$VCV[,'animal']/
  (heat_shrubt$VCV[,'animal']+heat_shrubt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_restricted <- fire_cues_heat3 %>%
  filter(Distribution == "Restricted")

n_distinct(heat_restricted$animal)
nrow(heat_restricted)

heat_restrictedt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                  random = ~ animal + ID + n:ID, 
                                  family = "multinomial2", pedigree = fire_cues_tree, 
                                  prior = priors, data = heat_restricted,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_restrictedt)
n_distinct(heat_restricted$animal)

lambda <- heat_restrictedt$VCV[,'animal']/
  (heat_restrictedt$VCV[,'animal']+heat_restrictedt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_widespread <- fire_cues_heat3 %>%
  filter(Distribution == "Widespread")

n_distinct(heat_widespread$animal)
nrow(heat_widespread)

heat_widespreadt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = fire_cues_tree, 
                                       prior = priors, data = heat_widespread,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_widespreadt)
n_distinct(heat_widespread$animal)

lambda <- heat_widespreadt$VCV[,'animal']/
  (heat_widespreadt$VCV[,'animal']+heat_widespreadt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_mesic_xeric <- fire_cues_heat3 %>%
  filter(Microhabitat == "Mesic/Xeric")

n_distinct(heat_mesic_xeric$animal)
nrow(heat_mesic_xeric)

heat_mesic_xerict <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = fire_cues_tree, 
                                       prior = priors, data = heat_mesic_xeric,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_mesic_xerict)
n_distinct(heat_mesic_xeric$animal)

lambda <- heat_mesic_xerict$VCV[,'animal']/
  (heat_mesic_xerict$VCV[,'animal']+heat_mesic_xerict$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_xeric <- fire_cues_heat3 %>%
  filter(Microhabitat == "Xeric")

n_distinct(heat_xeric$animal)
nrow(heat_xeric)

heat_xerict <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "multinomial2", pedigree = fire_cues_tree, 
                                        prior = priors, data = heat_xeric,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_xerict)
n_distinct(heat_xeric$animal)

lambda <- heat_xerict$VCV[,'animal']/
  (heat_xerict$VCV[,'animal']+heat_xerict$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_dormant <- fire_cues_heat3 %>%
  filter(Dormancy == "D")

n_distinct(heat_dormant$animal)
nrow(heat_dormant)

heat_dormantt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                  random = ~ animal + ID + n:ID, 
                                  family = "multinomial2", pedigree = fire_cues_tree, 
                                  prior = priors, data = heat_dormant,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_dormantt)
n_distinct(heat_dormant$animal)

lambda <- heat_dormantt$VCV[,'animal']/
  (heat_dormantt$VCV[,'animal']+heat_dormantt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_nondormant <- fire_cues_heat3 %>%
  filter(Dormancy == "ND")

n_distinct(heat_nondormant$animal)
nrow(heat_nondormant)

heat_nondormantt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                    random = ~ animal + ID + n:ID, 
                                    family = "multinomial2", pedigree = fire_cues_tree, 
                                    prior = priors, data = heat_nondormant,
                                    nitt = nite, thin = nthi, burnin = nbur, 
                                    verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                    pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_nondormantt)
n_distinct(heat_nondormant$animal)

lambda <- heat_nondormantt$VCV[,'animal']/
  (heat_nondormantt$VCV[,'animal']+heat_nondormantt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_mass <- fire_cues_heat3 %>%
  filter(!is.na(Dry_mass))

n_distinct(heat_mass$animal)
nrow(heat_mass)

heat_masst <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application * scale(Dry_mass), 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = fire_cues_tree, 
                                       prior = priors, data = heat_mass,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_masst)
n_distinct(heat_mass$animal)

lambda <- heat_masst$VCV[,'animal']/
  (heat_masst$VCV[,'animal']+heat_masst$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

#### Smoke ####

smoke_herbs <- fire_cues_smoke %>%
  filter(Growth_form == "Herb")

n_distinct(smoke_herbs$animal)
nrow(smoke_herbs)

smoke_herbst <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                  random = ~ animal + ID + n:ID, 
                                  family = "multinomial2", pedigree = fire_cues_tree, 
                                  prior = priors, data = smoke_herbs,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_herbst)
n_distinct(smoke_herbs$animal)

lambda <- smoke_herbst$VCV[,'animal']/
  (smoke_herbst$VCV[,'animal']+smoke_herbst$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

smoke_shrub <- fire_cues_smoke %>%
  filter(Growth_form == "Shrub")

n_distinct(smoke_shrub$animal)
nrow(smoke_shrub)

smoke_shrubt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                  random = ~ animal + ID + n:ID, 
                                  family = "multinomial2", pedigree = fire_cues_tree, 
                                  prior = priors, data = smoke_shrub,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_shrubt)
n_distinct(smoke_shrub$animal)

lambda <- smoke_shrubt$VCV[,'animal']/
  (smoke_shrubt$VCV[,'animal']+smoke_shrubt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

smoke_restricted <- fire_cues_smoke %>%
  filter(Distribution == "Restricted")

n_distinct(smoke_restricted$animal)
nrow(smoke_restricted)

smoke_restrictedt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = fire_cues_tree, 
                                       prior = priors, data = smoke_restricted,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_restrictedt)
n_distinct(smoke_restricted$animal)

lambda <- smoke_restrictedt$VCV[,'animal']/
  (smoke_restrictedt$VCV[,'animal']+smoke_restrictedt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

smoke_widespread <- fire_cues_smoke %>%
  filter(Distribution == "Widespread")

n_distinct(smoke_widespread$animal)
nrow(smoke_widespread)

smoke_widespreadt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = fire_cues_tree, 
                                       prior = priors, data = smoke_widespread,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_widespreadt)
n_distinct(smoke_widespread$animal)

lambda <- smoke_widespreadt$VCV[,'animal']/
  (smoke_widespreadt$VCV[,'animal']+smoke_widespreadt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

smoke_mesic_xeric <- fire_cues_smoke %>%
  filter(Microhabitat == "Mesic/Xeric")

n_distinct(smoke_mesic_xeric$animal)
nrow(smoke_mesic_xeric)

smoke_mesic_xerict <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                        random = ~ animal + n, 
                                        family = "multinomial2", pedigree = fire_cues_tree, 
                                        prior = priors2, data = smoke_mesic_xeric,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_mesic_xerict)
n_distinct(smoke_mesic_xeric$animal)

lambda <- smoke_mesic_xerict$VCV[,'animal']/
  (smoke_mesic_xerict$VCV[,'animal']+smoke_mesic_xerict$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

smoke_xeric <- fire_cues_smoke %>%
  filter(Microhabitat == "Xeric")

n_distinct(smoke_xeric$animal)
nrow(smoke_xeric)

smoke_xerict <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                  random = ~ animal + ID + n:ID, 
                                  family = "multinomial2", pedigree = fire_cues_tree, 
                                  prior = priors, data = smoke_xeric,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_xerict)
n_distinct(smoke_xeric$animal)

lambda <- smoke_xerict$VCV[,'animal']/
  (smoke_xerict$VCV[,'animal']+smoke_xerict$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

smoke_dormant <- fire_cues_smoke %>%
  filter(Dormancy == "D")

n_distinct(smoke_dormant$animal)
nrow(smoke_dormant)

smoke_dormantt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                    random = ~ animal + ID + n:ID, 
                                    family = "multinomial2", pedigree = fire_cues_tree, 
                                    prior = priors, data = smoke_dormant,
                                    nitt = nite, thin = nthi, burnin = nbur, 
                                    verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                    pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_dormantt)
n_distinct(smoke_dormant$animal)

lambda <- smoke_dormantt$VCV[,'animal']/
  (smoke_dormantt$VCV[,'animal']+smoke_dormantt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

smoke_nondormant <- fire_cues_smoke %>%
  filter(Dormancy == "ND")

n_distinct(smoke_nondormant$animal)
nrow(smoke_nondormant)

smoke_nondormantt <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "multinomial2", pedigree = fire_cues_tree, 
                                       prior = priors, data = smoke_nondormant,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_nondormantt)
n_distinct(smoke_nondormant$animal)

lambda <- smoke_nondormantt$VCV[,'animal']/
  (smoke_nondormantt$VCV[,'animal']+smoke_nondormantt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

smoke_mass <- fire_cues_smoke %>%
  filter(!is.na(Dry_mass))

n_distinct(smoke_mass$animal)
nrow(smoke_mass)

smoke_masst <- MCMCglmm::MCMCglmm(cbind(GermSeeds, nSeeds-GermSeeds) ~ Treatment_application * scale(Dry_mass), 
                                 random = ~ animal + ID + n:ID, 
                                 family = "multinomial2", pedigree = fire_cues_tree, 
                                 prior = priors, data = smoke_mass,
                                 nitt = nite, thin = nthi, burnin = nbur, 
                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                 pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_masst)
n_distinct(smoke_mass$animal)

lambda <- smoke_masst$VCV[,'animal']/
  (smoke_masst$VCV[,'animal']+smoke_masst$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

#### t50 ####

# Filtering the databse

heat_experiments_time <- germination %>%
  filter(ID %in% heat_studies, # Studies with whose IDs are in "heat_studies",
         Photoperiod != 0,
         !is.na(HeatShock_Temperature), # that they evaluated the effect of heat
         is.na(Smoke) | Smoke == "Control", # shocks alone (i.e., not in combination with smoke)
         !is.na(D1),
         !is.na(D2)) %>% 
  mutate(D0 = 0) %>%
  select(1:3, HeatShock_Temperature, HeatShock_Time, Smoke, GermSeeds, nSeeds, D0, D1:D30) %>%
  left_join(traits) %>%
  select(1:39, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits info.
  left_join(lcvp_names_checked) %>%
  filter(HeatShock_Temperature == "Control" | HeatShock_Temperature %in% c(80, 100, 200),
         Growth_form == "Herb" | Growth_form == "Shrub",
         Dormancy != "NC") %>%
  unite("HeatShock_Treatment", HeatShock_Temperature, HeatShock_Time, sep = "*", remove = FALSE) %>%
  mutate(HeatShock_Treatment = case_when(HeatShock_Treatment == "Control*Control" ~ "Control",
                                         HeatShock_Treatment != "Control*Control" ~ "Treatment"))

exclude_heat_experiments_time <- heat_experiments_time %>%
  select(1:6, GermSeeds, nSeeds) %>%
  mutate(GermProp = GermSeeds/nSeeds) %>%
  group_by(ID, Species_reported, Species_accepted, HeatShock_Treatment) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = HeatShock_Treatment, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 | Treatment < 0.1)) %>%
  select(1:3)

heat_experiments_time2 <- heat_experiments_time %>%
  anti_join(exclude_heat_experiments_time) %>% 
  unite("HeatShock_Treatment", HeatShock_Temperature, HeatShock_Time, sep = "*", remove = FALSE) %>%
  mutate(Treatment_type = "Single",
         Treatment = "Heat",
         HeatShock_Treatment = case_when(HeatShock_Treatment == "Control*Control" ~ "Control",
                                         HeatShock_Treatment != "Control*Control" ~ HeatShock_Treatment),
         Treatment_application = case_when(HeatShock_Treatment == "Control" ~ 0,
                                           .default = 1)) %>%
  filter(!(Species_acceptedLCVP == "Vellozia resinosa" & HeatShock_Treatment == "100*5")) %>%
  select(1:3, Treatment_type, Treatment, Treatment_application, HeatShock_Treatment, Smoke:Species_acceptedLCVP)

### Smoke

# Now, we are filtering the germination database: 

smoke_experiments_time <- germination %>%
  filter(ID %in% smoke_studies, # Studies with IDs in "smoke_studies",
         Photoperiod != 0,
         !is.na(D1),
         !is.na(D2),
         Smoke == "Smoke water" | Smoke == "Control",
         is.na(HeatShock_Temperature) | HeatShock_Temperature == "Control") %>%  # Assessed only the effect of smoke
  mutate(D0 = 0) %>%
  select(1:3, HeatShock_Temperature:Smoke, GermSeeds, nSeeds, D0, D1:D30) %>%
  left_join(traits) %>%
  select(1:44, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits info.
  left_join(lcvp_names_checked) %>%
  filter(Growth_form == "Herb" | Growth_form == "Shrub",
         Dormancy != "NC") %>%
  mutate(HeatShock_Treatment = NA) %>%
  select(1:3, HeatShock_Treatment, 4:39, Growth_form:Species_acceptedLCVP)

exclude_smoke_experiments_time <- smoke_experiments_time %>%
  select(1:3, Smoke, GermSeeds, nSeeds) %>%
  mutate(GermProp = GermSeeds/nSeeds,
         Smoke = case_when(Smoke == "Smoke water" ~ "Treatment",
                           .default = "Control")) %>%
  group_by(ID, Species_reported, Species_accepted, Smoke) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Smoke, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 | Treatment < 0.1)) %>%
  select(1:3)

smoke_experiments_time2 <- smoke_experiments_time %>%
  anti_join(exclude_smoke_experiments_time) %>%
  mutate(Treatment_type = "Single",
         Treatment = "Smoke",
         Treatment_application = case_when(Smoke == "Control" ~ 0,
                                           .default = 1)) %>%
  select(1:3, Treatment_type, Treatment, Treatment_application, HeatShock_Treatment, Smoke:Species_acceptedLCVP)

# Heat shock and smoke combined

heat_smoke_experiments_time <- germination %>%
  filter(ID %in% heat_smoke_studies,
         Photoperiod != 0,
         Smoke == "Smoke water" | Smoke == "Control",
         !is.na(HeatShock_Temperature),
         !is.na(D1),
         !is.na(D2)) %>%
  mutate(D0 = 0) %>%
  select(1:3, HeatShock_Temperature:Smoke, GermSeeds, nSeeds, D0, D1:D30) %>%
  left_join(traits) %>%
  select(1:45, Growth_form, Distribution, Microhabitat, Dry_mass, Dormancy) %>% # to recover the traits info.
  left_join(lcvp_names_checked) %>%
  filter(HeatShock_Temperature == "Control" | HeatShock_Temperature %in% c(80, 100, 200),
         Growth_form == "Herb" | Growth_form == "Shrub",
         Dormancy != "NC") %>%
  unite("Treatment", HeatShock_Temperature, HeatShock_Time, Smoke, sep = "*", remove = FALSE) %>%
  mutate(Treatment = case_when(Treatment == "Control*Control*Control" ~ "Control",
                               Treatment != "Control*Control*Control" ~ "Treatment")) %>%
  select(1:3, Treatment, 4:40, Growth_form:Species_acceptedLCVP)
  

exclude_heat_smoke_experiments_time <- heat_smoke_experiments_time %>%
  select(1:3, Treatment, GermSeeds, nSeeds) %>%
  mutate(GermProp = GermSeeds/nSeeds) %>%
  group_by(ID, Species_reported, Species_accepted, Treatment) %>%
  summarise(meanGermProp = mean(GermProp)) %>%
  tidyr::pivot_wider(names_from = Treatment, values_from = meanGermProp) %>%
  filter(is.na(Treatment) | is.na(Control) | (Control < 0.1 | Treatment < 0.1)) %>%
  select(1:3)

heat_smoke_experiments_time2 <- heat_smoke_experiments_time %>%
  anti_join(exclude_heat_smoke_experiments_time) %>%
  unite("HeatShock_Treatment", HeatShock_Temperature, HeatShock_Time, sep = "*", remove = FALSE) %>%
  mutate(Treatment_type = "Combined",
         Treatment = "Combined",
         Treatment_application = case_when(HeatShock_Treatment == "Control*Control" ~ 0,
                                           .default = 1)) %>%
  select(1:3, Treatment_type, Treatment, Treatment_application, HeatShock_Treatment, Smoke:Species_acceptedLCVP)

fire_cues_time_list <- list(heat_experiments_time2, smoke_experiments_time2, heat_smoke_experiments_time2)

fire_cues_time <- fire_cues_time_list %>%
  reduce(full_join) %>%
  filter(!is.na(Species_acceptedLCVP))

fire_cues_t50 <- germination.indices(as.data.frame(fire_cues_time),
                                     total.seeds.col = "nSeeds",
                                     counts.intervals.cols = intervals30,
                                     intervals = 1:31,
                                     t50 = TRUE,
                                     EmergenceRateIndex = FALSE,
                                     PeakGermPercent = FALSE,
                                     PeakGermTime = FALSE,
                                     GermValue = FALSE) %>%
  dplyr::select(1:47, t50_Farooq)
  

write.csv(fire_cues_t50, "fire_cues_t50_missing.csv")

fire_cues_time_final <- read_csv(str_c(git_find(),
                                       "/Analyses/Meta_Analyses/Fire/fire_cues_time_final.csv")) %>%
  mutate(ID = as.factor(ID),
         Treatment_application = as.factor(Treatment_application),
         n = row_number(),
         n = as.factor(n)) %>%
  rename(animal = Species_acceptedLCVP)

n_distinct(fire_cues_time_final$animal)
nrow(fire_cues_time_final)

fire_cues_time_names <- fire_cues_time_final %>% 
  dplyr::select(animal) %>%
  distinct() %>%
  column_to_rownames(var = "animal")

fire_cues_time_names <- names_for_pruning(fire_cues_time_names)

# Now, we prune the phylogeny in order to keep only these names.

fire_cues_time_tree <- prune.sample(fire_cues_time_names, phylo = full_tree)
fire_cues_time_tree$node.label <- NULL

fire_cues_time_all <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                    random = ~ animal + ID + n:ID, 
                                    family = "gaussian", pedigree = fire_cues_time_tree, 
                                    prior = priors, data = fire_cues_time_final,
                                    nitt = nite, thin = nthi, burnin = nbur, 
                                    verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                    pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_all)
n_distinct(fire_cues_time_final$animal)

lambda <- fire_cues_time_all$VCV[,'animal']/
  (fire_cues_time_all$VCV[,'animal']+fire_cues_time_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_time_combined <- fire_cues_time_final %>%
  filter(Treatment_type == "Combined")

n_distinct(fire_cues_time_combined$animal)
nrow(fire_cues_time_combined)

fire_cues_time_combinedt <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                          random = ~ animal + n, 
                                          family = "gaussian", pedigree = fire_cues_time_tree, 
                                          prior = priors2, data = fire_cues_time_combined,
                                          nitt = nite, thin = nthi, burnin = nbur, 
                                          verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                          pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_combinedt)
n_distinct(fire_cues_time_combinedt$animal)

lambda <- fire_cues_time_combinedt$VCV[,'animal']/
  (fire_cues_time_combinedt$VCV[,'animal']+fire_cues_time_combinedt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_time_single <- fire_cues_time_final %>%
  filter(Treatment_type == "Single")

n_distinct(fire_cues_time_single$animal)
nrow(fire_cues_time_single)

fire_cues_time_singlet <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "gaussian", pedigree = fire_cues_time_tree, 
                                        prior = priors, data = fire_cues_time_single,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_singlet)
n_distinct(fire_cues_time_single$animal)

lambda <- fire_cues_time_singlet$VCV[,'animal']/
  (fire_cues_time_singlet$VCV[,'animal']+fire_cues_time_singlet$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_time_heat <- fire_cues_time_final %>%
  filter(Treatment == "Heat") %>%
  mutate(HeatShock_Treatment = as.factor(HeatShock_Treatment))

n_distinct(fire_cues_time_heat$animal)
nrow(fire_cues_time_heat)

fire_cues_time_heatt <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                      random = ~ animal + ID + n:ID, 
                                      family = "gaussian", pedigree = fire_cues_time_tree, 
                                      prior = priors, data = fire_cues_time_heat,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_heatt)
n_distinct(fire_cues_time_final$animal)

lambda <- fire_cues_time_heatt$VCV[,'animal']/
  (fire_cues_time_heatt$VCV[,'animal']+fire_cues_time_heatt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_time_smoke <- fire_cues_time_final %>%
  filter(Treatment == "Smoke")

n_distinct(fire_cues_time_smoke$animal)
nrow(fire_cues_time_smoke)

fire_cues_time_smoket <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "gaussian", pedigree = fire_cues_time_tree, 
                                       prior = priors, data = fire_cues_time_smoke,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_smoket)
n_distinct(fire_cues_time_smoke$animal)

lambda <- fire_cues_time_smoket$VCV[,'animal']/
  (fire_cues_time_smoket$VCV[,'animal']+fire_cues_time_smoket$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_time_dormant <- fire_cues_time_final %>%
  filter(Dormancy == "D")

n_distinct(fire_cues_time_dormant$animal)
nrow(fire_cues_time_dormant)

fire_cues_time_dormantt <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                         random = ~ animal + ID + n:ID, 
                                         family = "gaussian", pedigree = fire_cues_time_tree,
                                         prior = priors, data = fire_cues_time_dormant,
                                         nitt = nite, thin = nthi, burnin = nbur, 
                                         verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                         pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_dormantt)
n_distinct(fire_cues_time_dormant$animal)

lambda <- fire_cues_time_dormantt$VCV[,'animal']/
  (fire_cues_time_dormantt$VCV[,'animal']+fire_cues_time_dormantt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_time_nondormant <- fire_cues_time_final %>%
  filter(Dormancy == "ND")

n_distinct(fire_cues_time_nondormant$animal)
nrow(fire_cues_time_nondormant)

fire_cues_time_non_dormantt <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                                  random = ~ animal + ID + n:ID, 
                                                  family = "gaussian", pedigree = fire_cues_time_tree,
                                             prior = priors, data = fire_cues_time_nondormant,
                                             nitt = nite, thin = nthi, burnin = nbur, 
                                             verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                             pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_non_dormantt)
n_distinct(fire_cues_nondormant$animal)

lambda <- fire_cues_time_non_dormantt$VCV[,'animal']/
  (fire_cues_time_non_dormantt$VCV[,'animal']+fire_cues_time_non_dormantt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_time_mass <- fire_cues_time_final %>%
  filter(!is.na(Dry_mass))

n_distinct(fire_cues_time_mass$animal)
nrow(fire_cues_time_mass)

fire_cues_time_masst <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application * scale(Dry_mass), 
                                           random = ~ animal + ID + n:ID, 
                                           family = "gaussian", pedigree = fire_cues_time_tree,
                                      prior = priors, data = fire_cues_time_mass,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_masst)
n_distinct(fire_cues_mass$animal)

lambda <- fire_cues_time_masst$VCV[,'animal']/
  (fire_cues_time_masst$VCV[,'animal']+fire_cues_time_masst$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Now we will test the effect of different heat shock temperatures

fire_cues_time_heat2 <- fire_cues_time_heat %>%
  filter(HeatShock_Treatment != "100*3") # We are excluding this treatment as it was only tested in two species in a      single study.

fire_cues_time_80_2 <- fire_cues_time_heat2 %>%
  filter(HeatShock_Treatment == "Control" | HeatShock_Treatment == "80*2")

n_distinct(fire_cues_time_80_2$animal)
nrow(fire_cues_time_80_2)

fire_cues_time_80_2t <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                           random = ~ animal + n, 
                                           family = "gaussian", pedigree = fire_cues_time_tree,
                                      prior = priors2, data = fire_cues_time_80_2,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_80_2t)
n_distinct(figure_cues_80_2$animal)

lambda <- fire_cues_time_80_2t$VCV[,'animal']/
  (fire_cues_time_80_2t$VCV[,'animal']+fire_cues_time_80_2t$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

fire_cues_time_80_5 <- fire_cues_time_heat2 %>%
  filter(HeatShock_Treatment == "Control" | HeatShock_Treatment == "80*5")

n_distinct(fire_cues_time_80_5$animal)
nrow(fire_cues_time_80_5)

fire_cues_time_80_5t <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                           random = ~ animal + n, 
                                           family = "gaussian", pedigree = fire_cues_time_tree,
                                           prior = priors2, data = fire_cues_time_80_5,
                                      nitt = nite, thin = nthi, burnin = nbur, 
                                      verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                      pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_80_5t)
n_distinct(figure_cues_80_5$animal)

lambda <- fire_cues_time_80_5t$VCV[,'animal']/
  (fire_cues_time_80_5t$VCV[,'animal']+fire_cues_time_80_5t$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)


fire_cues_time_100_5 <- fire_cues_time_heat2 %>%
  filter(HeatShock_Treatment == "Control" | HeatShock_Treatment == "100*5")

n_distinct(fire_cues_time_100_5$animal)
nrow(fire_cues_time_100_5)

fire_cues_time_100_5t <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                            random = ~ animal + n, 
                                            family = "gaussian", pedigree = fire_cues_time_tree,
                                            prior = priors2, data = fire_cues_time_100_5,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(fire_cues_time_100_5t)
n_distinct(figure_cues_100_5$animal)

lambda <- fire_cues_time_100_5t$VCV[,'animal']/
  (fire_cues_time_100_5t$VCV[,'animal']+fire_cues_time_100_5t$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

# Now, we'll be testing the effect of heat shocks on the different plant groups

fire_cues_time_heat3 <- fire_cues_time_heat %>%
  filter(HeatShock_Treatment != "200*1")

n_distinct(fire_cues_time_heat3$animal)
nrow(fire_cues_time_heat3)

heat_time_all <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                               random = ~ animal + ID + n:ID, 
                               family = "gaussian", pedigree = fire_cues_time_tree,
                               prior = priors, data = fire_cues_time_heat3,
                               nitt = nite, thin = nthi, burnin = nbur, 
                               verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                               pr = FALSE, pl = FALSE)

beepr::beep()


summary(heat_time_all)
n_distinct(fire_cues_time_heat3$animal)

lambda <- heat_time_all$VCV[,'animal']/
  (heat_time_all$VCV[,'animal']+heat_time_all$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_time_herbs <- fire_cues_time_heat3 %>%
  filter(Growth_form == "Herb")

n_distinct(heat_time_herbs$animal)
nrow(heat_time_herbs)

heat_time_herbst <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "gaussian", pedigree = fire_cues_time_tree,
                                  prior = priors, data = heat_time_herbs,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_time_herbst)
n_distinct(heat_time_herbs$animal)

lambda <- heat_time_herbst$VCV[,'animal']/
  (heat_time_herbst$VCV[,'animal']+heat_time_herbst$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_time_shrub <- fire_cues_time_heat3 %>%
  filter(Growth_form == "Shrub")

n_distinct(heat_time_shrub$animal)
nrow(heat_time_shrub)

heat_time_shrubt <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "gaussian", pedigree = fire_cues_time_tree,
                                  prior = priors, data = heat_time_shrub,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_time_shrubt)
n_distinct(heat_shrub$animal)

lambda <- heat_time_shrubt$VCV[,'animal']/
  (heat_time_shrubt$VCV[,'animal']+heat_time_shrubt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_time_restricted <- fire_cues_time_heat3 %>%
  filter(Distribution == "Restricted")

n_distinct(heat_time_restricted$animal)
nrow(heat_time_restricted)

heat_time_restrictedt <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = fire_cues_time_tree,
                                       prior = priors, data = heat_time_restricted,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_time_restrictedt)
n_distinct(heat_time_restricted$animal)

lambda <- heat_time_restrictedt$VCV[,'animal']/
  (heat_time_restrictedt$VCV[,'animal']+heat_time_restrictedt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_time_widespread <- fire_cues_time_heat3 %>%
  filter(Distribution == "Widespread")

n_distinct(heat_time_widespread$animal)
nrow(heat_time_widespread)

heat_time_widespreadt <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = fire_cues_time_tree,
                                       prior = priors, data = heat_time_widespread,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_time_widespreadt)
n_distinct(heat_widespread$animal)

lambda <- heat_time_widespreadt$VCV[,'animal']/
  (heat_time_widespreadt$VCV[,'animal']+heat_time_widespreadt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_time_mesic_xeric <- fire_cues_time_heat3 %>%
  filter(Microhabitat == "Mesic/Xeric")

n_distinct(heat_time_mesic_xeric$animal)
nrow(heat_time_mesic_xeric)

heat_time_mesic_xerict <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                             random = ~ animal + ID + n:ID, 
                                             family = "gaussian", pedigree = fire_cues_time_tree,
                                        prior = priors, data = heat_time_mesic_xeric,
                                        nitt = nite, thin = nthi, burnin = nbur, 
                                        verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                        pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_time_mesic_xerict)
n_distinct(heat_time_mesic_xeric$animal)

lambda <- heat_time_mesic_xerict$VCV[,'animal']/
  (heat_time_mesic_xerict$VCV[,'animal']+heat_time_mesic_xerict$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_time_xeric <- fire_cues_time_heat3 %>%
  filter(Microhabitat == "Xeric")

n_distinct(heat_time_xeric$animal)
nrow(heat_time_xeric)

heat_time_xerict <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                       random = ~ animal + ID + n:ID, 
                                       family = "gaussian", pedigree = fire_cues_time_tree,
                                  prior = priors, data = heat_time_xeric,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_time_xerict)
n_distinct(heat_xeric$animal)

lambda <- heat_time_xerict$VCV[,'animal']/
  (heat_time_xerict$VCV[,'animal']+heat_time_xerict$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_time_dormant <- fire_cues_time_heat3 %>%
  filter(Dormancy == "D")

n_distinct(heat_time_dormant$animal)
nrow(heat_time_dormant)

heat_time_dormantt <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                         random = ~ animal + n, 
                                         family = "gaussian", pedigree = fire_cues_time_tree,
                                    prior = priors2, data = heat_time_dormant,
                                    nitt = nite, thin = nthi, burnin = nbur, 
                                    verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                    pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_time_dormantt)
n_distinct(heat_time_dormant$animal)

lambda <- heat_time_dormantt$VCV[,'animal']/
  (heat_time_dormantt$VCV[,'animal']+heat_time_dormantt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_time_nondormant <- fire_cues_time_heat3 %>%
  filter(Dormancy == "ND")

n_distinct(heat_time_nondormant$animal)
nrow(heat_time_nondormant)

heat_time_nondormantt <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                            random = ~ animal + ID + n:ID, 
                                            family = "gaussian", pedigree = fire_cues_time_tree,
                                       prior = priors, data = heat_time_nondormant,
                                       nitt = nite, thin = nthi, burnin = nbur, 
                                       verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                       pr = FALSE, pl = FALSE)

beepr::beep()

summary(heat_nondormantt)
n_distinct(heat_nondormant$animal)

lambda <- heat_time_nondormantt$VCV[,'animal']/
  (heat_time_nondormantt$VCV[,'animal']+heat_time_nondormantt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

heat_time_mass <- fire_cues_time_heat3 %>%
  filter(!is.na(Dry_mass))

n_distinct(heat_time_mass$animal)
nrow(heat_time_mass)

heat_time_masst <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application * scale(Dry_mass), 
                                      random = ~ animal + ID + n:ID, 
                                      family = "gaussian", pedigree = fire_cues_time_tree,
                                 prior = priors, data = heat_time_mass,
                                 nitt = nite, thin = nthi, burnin = nbur, 
                                 verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                 pr = FALSE, pl = FALSE)
beepr::beep()

summary(heat_time_masst)
n_distinct(heat_time_mass$animal)

lambda <- heat_time_masst$VCV[,'animal']/
  (heat_time_masst$VCV[,'animal']+heat_time_masst$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

#### Smoke ####

smoke_time_herbs <- fire_cues_time_smoke %>%
  filter(Growth_form == "Herb")

n_distinct(smoke_time_herbs$animal)
nrow(smoke_time_herbs)

smoke_time_herbst <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                        random = ~ animal + ID + n:ID, 
                                        family = "gaussian", pedigree = fire_cues_time_tree,
                                   prior = priors, data = smoke_time_herbs,
                                   nitt = nite, thin = nthi, burnin = nbur, 
                                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                   pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_time_herbst)
n_distinct(smoke_herbs$animal)

lambda <- smoke_time_herbst$VCV[,'animal']/
  (smoke_time_herbst$VCV[,'animal']+smoke_time_herbst$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

smoke_time_shrub <- fire_cues_time_smoke %>%
  filter(Growth_form == "Shrub")

n_distinct(smoke_time_shrub$animal)
nrow(smoke_time_shrub)

smoke_time_shrubt <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application, 
                                   random = ~ animal + ID + n:ID, 
                                   family = "gaussian", pedigree = fire_cues_time_tree,
                                   prior = priors, data = smoke_time_shrub,
                                   nitt = nite, thin = nthi, burnin = nbur, 
                                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                   pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_time_shrubt)
n_distinct(smoke_shrub$animal)

lambda <- smoke_time_shrubt$VCV[,'animal']/
  (smoke_time_shrubt$VCV[,'animal']+smoke_time_shrubt$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

smoke_time_mass <- fire_cues_time_smoke %>%
  filter(!is.na(Dry_mass))

n_distinct(smoke_time_mass$animal)
nrow(smoke_time_mass)

smoke_time_masst <- MCMCglmm::MCMCglmm(t50_Farooq ~ Treatment_application * scale(Dry_mass), 
                                  random = ~ animal + ID + n:ID, 
                                  family = "gaussian", pedigree = fire_cues_time_tree, 
                                  prior = priors, data = smoke_time_mass,
                                  nitt = nite, thin = nthi, burnin = nbur, 
                                  verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                                  pr = FALSE, pl = FALSE)

beepr::beep()

summary(smoke_time_masst)
n_distinct(smoke_mass$animal)

lambda <- smoke_time_masst$VCV[,'animal']/
  (smoke_time_masst$VCV[,'animal']+smoke_time_masst$VCV[,'units'])

posterior.mode(lambda)
HPDinterval(lambda)

