# The following code fits phylogenetic generalised least squares (PGLS) and phylogenetic logistic regression (PLR) models to test the variation of seed traits between different ecological groups (i.e., growth forms, microhabitat, geographic distributions).

# Author: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com)
# Last update: February 29th, 2024.

#### 0. Loading packages (assuming they are already installed) ####

library(gert) # v. 2.0.1 for managing Git repositories.
library(stringr) # v. 1.5.1 for string manipulation.
library(readr) # v. 2.1.5 for loading .csv files.
library(dplyr) # v. 1.1.4 for database manipulation.
library(tibble) # v. 3.2.1 for dataframe manipulations.
library(caper) # v. 1.0.3 for PGLS implementation.
library(picante) # v. 1.8.2. for pruning our phylogeny.
library(car) # v. 3.1-2 for logit transformation.
library(phylolm) # v. 2.6.2 for phylogenetic logistic regression.

# We are also loading the R-package daee (v. 0.1.7), available on GitHub.

library(devtools) # v. 2.4.5 for downloading packages from GitHub.

install_github("vanderleidebastiani/daee")
library(daee)

#### 1. Set working directory ####

setwd(str_c(git_find(),"/Analyses/Trait_Variation"))

#### 2. Loading files ####

# First, the phylogeny we created (see Analyses/Phylogeny/phylogeny.R).

full_tree <- read.tree(str_c(git_find(),"/Analyses/Phylogeny/outcrop_phylo.tre"))

# Second, we will upload the "traits_joined.csv" dataset that we created while evaluating the phylogenetic signal (see Analyses/Phylogenetic_Signal).

read_csv((str_c(git_find(),"/Analyses/Phylogenetic_Signal/traits_joined.csv")))

#### 3. Phylogenetic generalised least squares (PGLS) models for seed mass ####

# Now that we've assessed the phylogenetic signal of our functional traits, we will test if these traits vary between the ecological groups we have defined according to species growth form, geographic distribution and microhabitat.

# For this purpose, we are using a pgls as implemented in the R package caper. We need two objects: i) the phylogeny (which we already created while calculating the phylogenetic signal of our traits) and ii) an special object called "comparative.data" object.

# First, we are creating a dataframe with the mean seed mass of each species, and its growth form, distribution and microhabitat.

mass_traits <- traits_joined %>%
  dplyr::select(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat, Dry_mass) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Dry_mass),
         Growth_form != "Liana",
         Growth_form != "Succulent",
         Growth_form != "Tree") %>%
  group_by(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat) %>%
  summarise(Mass = mean(Dry_mass))

# Preliminary checks in the database showed that few species (6) had records for more than one growth form, microhabitat or geographic distribution class. In these cases we recorded microhabitat as 'Mesic/Xeric' and the growth form reported in Flora do Brasil.These changes were made my manually modifying the following database.

write.csv(mass_traits, file = "check_mass_traits.csv", row.names = FALSE)

# After doing these changes, we are uploading again the dataset.

mass_traits <- read.csv(str_c(git_find(),"/Analyses/Trait_Variation/checked_mass_traits.csv"))

mass_traits2 <- mass_traits %>%
  column_to_rownames(var = "Species_acceptedLCVP")

# Now, we are taking the names in the mass_traits database to prune the phylogenetic tree. This procedure is the same as carried for the phylogenetic signal section. For that purpose we will use the names_for_pruning function we created in Phylogenetic_Signal.R

mass_traits_names <- names_for_pruning(mass_traits2)

# Now that we have the object with the names, we'll prune the full phylogenetic tree to keep only the species that have seed mass data.

mass_traits_tree <- prune.sample(mass_traits_names, phylo = full_tree)
mass_traits_tree$node.label <- NULL

# Wtih this prune phylogeny, we'll create the comparative.data object, using the in built function comparative.data in the caper package.

mass_comparative <- comparative.data(phy = mass_traits_tree,
                                     data = mass_traits,
                                     names.col = 'Species_acceptedLCVP',
                                     vcv = TRUE,
                                     na.omit = FALSE)

# With these two objects created, we proceed to run the pgls, by first analyzing a simple additive model

mass_pgls <- pgls(log(Mass) ~ Growth_form + Distribution + Microhabitat,
                 data = mass_comparative, 
                 lambda = "ML")

# Before seeing the results of the model, we are checking the normality of the residuals and whether these are non-independent.

mass_pgls_residuals <- data.frame(mass_pgls$residuals)

ks.test(mass_pgls_residuals$mass_pgls.residuals, pnorm,
        mean(mass_pgls_residuals$mass_pgls.residuals), 
        sd(mass_pgls_residuals$mass_pgls.residuals))

# Residuals are normally distributed.

plot(mass_pgls)

# Fitted values seem independent.

# Now we proceed to see the results of the model.

summary(mass_pgls)

# The model suggest that the data provide evidence of seed mass variation between growth forms and microhabitats, with shrubs and species from Mesic/Xeric habitats producing heavier seeds. Now, we are testing a model considering all possible interaction between predictors.

mass_pgls2 <- pgls(log(Mass) ~ Growth_form * Distribution * Microhabitat,
                   data = mass_comparative,
                   lambda = "ML")

# We are testing the normality of the residuals...

mass_pgls_residuals2 <- data.frame(mass_pgls2$residuals)

ks.test(mass_pgls_residuals2$mass_pgls2.residuals, pnorm,
        mean(mass_pgls_residuals2$mass_pgls2.residuals), 
        sd(mass_pgls_residuals2$mass_pgls2.residuals))

# and visually inspecting the distribution of fitted values. Both are OK.

plot(mass_pgls2)

# Now we are checking the signifance of each predictor.

anova.pgls(mass_pgls2)

# Since none of the interactions between predictors were statistically significant, our analysis were restricted to the additive model. Moreover, the Akaike Information Criteria (AIC) suggest that the simple additive model performs slightly better.

AIC(mass_pgls, mass_pgls2)

#### 4. PGLS for seed water content ####

# Now, we repeat the same routine with the seed water content.

water_traits <- traits_joined %>%
  dplyr::select(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat, Water_content) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Water_content),
         Growth_form != "Liana",
         Growth_form != "Succulent",
         Growth_form != "Tree") %>%
  group_by(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat) %>%
  summarise(Water = mean(Water_content/100)) # We are transforming percentages into proportions so we can logit-transform the values latter.

# Unlike the seed mass database, we did not found any species with records for more than one class for any of the predictors. As a result, no manual changes were required.

water_traits <- data.frame(water_traits)

# We proceed to prune the phylogenetic tree to keep only the species in the water_traits database.

water_traits2 <- water_traits %>%
  column_to_rownames(var = "Species_acceptedLCVP")

water_traits_names <- names_for_pruning(water_traits2)

water_traits_tree <- prune.sample(water_traits_names, phylo = full_tree)
water_traits_tree$node.label <- NULL

# After pruning the tree, we create the comparative.data object.

water_comparative <- comparative.data(phy = water_traits_tree,
                                      data = water_traits,
                                      names.col = 'Species_acceptedLCVP',
                                      vcv = TRUE,
                                      na.omit = FALSE)

# Having these two objects, we run the additive model.

water_pgls <- pgls(logit(Water) ~ Growth_form + Distribution + Microhabitat,
                   data = water_comparative,
                   lambda = "ML")

# We check the assumptions of the model. Both seem OK.

water_pgls_residuals <- data.frame(water_pgls$residuals)

ks.test(water_pgls_residuals$water_pgls.residuals, pnorm,
        mean(water_pgls_residuals$water_pgls.residuals), 
        sd(water_pgls_residuals$water_pgls.residuals))

plot(water_pgls)

# Now we check the statistical significance of the predictor to the model.

summary(water_pgls)

# The model indicates that the data do not provide evidence of variation in seed water content between any of the tested ecological groups.

# Now we run the model, this time considering the potential interaction between predictors.

water_pgls2 <- pgls(logit(Water) ~ Growth_form * Distribution * Microhabitat,
                    data = water_comparative,
                    lambda = "ML")

# We check the assumptions of the model. Both seem OK.

water_pgls_residuals2 <- data.frame(water_pgls2$residuals)

ks.test(water_pgls_residuals2$water_pgls2.residuals, pnorm,
        mean(water_pgls_residuals2$water_pgls2.residuals), 
        sd(water_pgls_residuals2$water_pgls2.residuals))

plot(water_pgls2)

# Now we check the statistical significance of the predictor to the model.

summary(water_pgls2)

# Since any of the intractions between predictors were stasitical significant, we restricted our analysis to the additive model. Also, the AIC indicates that the simple additive models is better.

AIC(water_pgls, water_pgls2)

#### 5. PGLS for percentage of embryoless seeds ####

# We are repeating the same steps with the percentage of embryoless seeds.

embryoless_traits <- traits_joined %>%
  dplyr::select(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat, Embryoless) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Embryoless),
         Growth_form != "Liana",
         Growth_form != "Succulent",
         Growth_form != "Tree") %>%
  group_by(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat) %>%
  summarise(Embryoless = mean(Embryoless/100)) # We are transforming percentages intro proportions.

# Unlike the seed mass database, we did not found any species with records for more than one class for any of the predictors. As a result, no manual changes were required.

# We proceed with pruning the phylogenetic tree and creating the comparative.data object.

embryoless_traits <- data.frame(embryoless_traits)

embryoless_traits2 <- embryoless_traits %>%
  column_to_rownames(var = "Species_acceptedLCVP")

embryoless_traits_names <- names_for_pruning(embryoless_traits2)
embryoless_traits_tree <- prune.sample(embryoless_traits_names, phylo = full_tree)
embryoless_traits_tree$node.label <- NULL

# Now we create the comparative.data object.

embryoless_comparative <- comparative.data(phy = embryoless_traits_tree,
                                           data = embryoless_traits,
                                           names.col = 'Species_acceptedLCVP',
                                           vcv = TRUE,
                                           na.omit = FALSE)

# Now, we proceed to run the additive model.

embryoless_pgls <- pgls(logit(Embryoless) ~ Growth_form + Distribution + Microhabitat,
                        data = embryoless_comparative,
                        lambda = "ML")

# We check the model assumptions of residuals normality and fitted values. Both are met.

embryoless_pgls_residuals <- data.frame(embryoless_pgls$residuals)

ks.test(embryoless_pgls_residuals$embryoless_pgls.residuals, pnorm,
        mean(embryoless_pgls_residuals$embryoless_pgls.residuals), 
        sd(embryoless_pgls_residuals$embryoless_pgls.residuals))

plot(embryoless_pgls)

# We check the statistical significance of  each predictor.

summary(embryoless_pgls)

# The model indicates that there is only weak evidence of variation on the percentage of embryoless seeds, with shrubs exhibiting more embryoless seeds than herbs.

# We re run the model considering the potential interactions between predictors.

embryoless_pgls2 <- pgls(logit(Embryoless) ~ Growth_form * Distribution * Microhabitat,
                         data = embryoless_comparative,
                         lambda = "ML")

# We check the model assumptions of residuals normality and fitted values. Both are met.

embryoless_pgls_residuals2 <- data.frame(embryoless_pgls2$residuals)

ks.test(embryoless_pgls_residuals2$embryoless_pgls2.residuals, pnorm,
        mean(embryoless_pgls_residuals2$embryoless_pgls2.residuals), 
        sd(embryoless_pgls_residuals2$embryoless_pgls2.residuals))

plot(embryoless_pgls2)

# We check the statistical significance of  each predictor.

anova(embryoless_pgls2)

# Since any of the interactions between predictors were statistical significant,we restricted our analysis to the additive model. Also, the AIC suggest that the simple additive model performs better.

AIC(embryoless_pgls, embryoless_pgls2)

#### 6. PGLS for percentage of viable seeds ####

# To conclude the analysis of quantiative seed traits, we repeat the routine for the percentage of embryoless seeds (%).

viable_traits <- traits_joined %>%
  dplyr::select(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat, Viable) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Viable),
         Growth_form != "Liana",
         Growth_form != "Succulent",
         Growth_form != "Tree") %>%
  group_by(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat) %>%
  summarise(Viable = mean(Viable/100)) # We convert the percentage into proportions to allow the logit transformation.

# Unlike the seed mass database, we did not found any species with records for more than one class for any of the predictors. As a result, no manual changes were required.

# Second, we prune the phylogenetic tree and the comparative.data object.

viable_traits <- data.frame(viable_traits)

viable_traits2 <- viable_traits %>%
  column_to_rownames(var = "Species_acceptedLCVP")

viable_traits_names <- names_for_pruning(viable_traits2)

viable_traits_tree <- prune.sample(viable_traits_names, phylo = full_tree)
viable_traits_tree$node.label <- NULL

viable_comparative <- comparative.data(phy = viable_traits_tree,
                                       data = viable_traits,
                                       names.col = 'Species_acceptedLCVP',
                                       vcv = TRUE,
                                       na.omit = FALSE)

# We proceed to run the additive model.

viable_pgls <- pgls(logit(Viable) ~ Growth_form + Distribution + Microhabitat,
                    data = viable_comparative,
                    lambda = "ML")

# Checking model assumptions. Both are met.

viable_pgls_residuals <- data.frame(viable_pgls$residuals)

ks.test(viable_pgls_residuals$viable_pgls.residuals, pnorm,
        mean(viable_pgls_residuals$viable_pgls.residuals), 
        sd(viable_pgls_residuals$viable_pgls.residuals))

plot(viable_pgls)

# We test for the statistical significance of predictors.

summary(viable_pgls)

# The model suggest that the data do not provide evidence of variation on the percentage of viable seeds across ecological groups.

# We re-run the model considering potential interactions between predictors.

viable_pgls2 <- pgls(logit(Viable) ~ Growth_form * Distribution * Microhabitat,
                     data = viable_comparative,
                     lambda = "ML")

# Checking model assumptions. Since residuals were only marginally different from normality we assumed both assumptions were met.

viable_pgls_residuals2 <- data.frame(viable_pgls2$residuals)

ks.test(viable_pgls_residuals2$viable_pgls2.residuals, pnorm,
        mean(viable_pgls_residuals2$viable_pgls2.residuals), 
        sd(viable_pgls_residuals2$viable_pgls2.residuals))

plot(viable_pgls)

# We test for the statistical significance of predictors.

anova(viable_pgls2)

# Since none of the interactions hold stasitical significance, we restricted our analysis to the additive model. Moreover, the AIC of the simple additive model indicactes that it permofed better.

AIC(viable_pgls, viable_pgls2)

#### 7. Phylogenetic logistic regression (PLR) of seed dormancy. ####

# After assesing potential differences in quantitative traits, now we proceed to to doing the same analysis with qualitative traits, starting with seed dormancy. For that purpose, we implemented a phylogenetic logistic regression as implemented in the phylolm package.

# First, we create a dataset with the seed dormancy information and the traits we are interested in.

dormancy_traits <- traits_joined %>%
  dplyr::select(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat, Dormancy) %>%
  filter(!is.na(Species_acceptedLCVP),
         Dormancy != "NC",
         Growth_form != "Liana",
         Growth_form != "Succulent",
         Growth_form != "Tree") %>%
  group_by(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat, Dormancy) %>%
  distinct()

write.csv(dormancy_traits, "check_dormancy_traits.csv")

# A preliminary check of the database showed that several species had reports of more than one class for at least one of the predictors or the dormancy variable. Twelve species that had records of producing both dormant and non-dormant seeds, each of these records was considered as different populations. Species with more than one growth-form reported were assigned the growth form reported at the Flora e Funga do Brasil site. This was the case of three species (Chamaecrista rotundifolia, Vellozia alutacea and V. epidendroides). Six species had records 
# in more than one microhabitat class. All these were assigned to the 'Mesic/Xeric' category.

dormancy_traits <- read.csv(str_c(git_find(),"/Analyses/Trait_Variation/checked_dormancy_traits.csv"))

dormancy_traits <- dormancy_traits %>%
  filter(!is.na(Growth_form), 
         !is.na(Distribution), 
         !is.na(Microhabitat)) %>%
  column_to_rownames("Species_acceptedLCVP")

# We recode the dormancy trait, with non-dormant (ND) seeds = 0, and dormant
# (D) = 1...

dormancy_traits$Dormancy <- dplyr::recode(dormancy_traits$Dormancy, "ND" = 0, "D" = 1)

# ... and create the object with the names for pruning the phylogeny.

dormancy_traits_names <- names_for_pruning(dormancy_traits)

# Now, we load a dataset that is going to allow us to bind our "populations" to their respective species. The dataset contains three columns. The first is the species which the population is going to be binded, the second is the population to be binded. The third column in the length of the branch we are adding. For this case, we are leaving it as NA, so the function uses a branch length that keeps the tree ultrametric.

populations_dormancy_traits <- read.csv(str_c(git_find(),
                                              "/Analyses/Trait_Variation/populations_dormancy_traits.csv"))
populations_dormancy_traits <- as.matrix(populations_dormancy_traits)
colnames(populations_dormancy_traits) <- NULL
populations_dormancy_traits

# We are adding this populations to the full phylogenetic tree...

dormancy_traits_tree <- add.taxa.phylo(full_tree, populations_dormancy_traits)

# and pruning it, so that it contains only species with seed dormancy data.

dormancy_traits_tree2 <- prune.sample(dormancy_traits_names, phylo = dormancy_traits_tree$tree)

# We are going to confirm that the tree tips and the species in the dataset are in the same order.

cbind(row.names(dormancy_traits), dormancy_traits_tree2$tip.label)
dormancy_traits <- dormancy_traits[order(match(rownames(dormancy_traits), 
                                               dormancy_traits_tree2$tip.label)), , drop = FALSE]
cbind(row.names(dormancy_traits), dormancy_traits_tree2$tip.label)

# With all these objects created, we are running the phylogenetic logistic regression. We start with a simple additive model.

dormancy_plr <- phyloglm(Dormancy ~ Growth_form + Distribution + Microhabitat,
                         data = dormancy_traits, phy = dormancy_traits_tree2)

summary(dormancy_plr)

# And we also test a model considering all possible interactions between traits.

dormancy_plr2 <- phyloglm(Dormancy ~ Growth_form * Distribution * Microhabitat,
                          data = dormancy_traits, phy = dormancy_traits_tree2)

summary(dormancy_plr2)

# Since the interactions between predictors did not yield a significant estimate, we restricted our analysis to the simple additive model.Also, the AIC of the simple additive model indicate that they performed slightly better.

AIC(dormancy_plr)
AIC(dormancy_plr2)

#### 8. PLR for dispersal syndromes ####

# Now we repeat the same routine with seed dispersal syndrome.

# We build a dataset with the seed dispersal syndrome information.

dispersal_traits <- traits_joined %>%
  dplyr::select(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat, Dispersal_syndrome) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Dispersal_syndrome),
         Growth_form != "Liana",
         Growth_form != "Succulent",
         Growth_form != "Tree") %>%
  group_by(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat, Dispersal_syndrome) %>%
  distinct()

write.csv(dispersal_traits, "check_dispersal_traits.csv")

# Now, we are looking on whether there are multiple combinations of species traits and predictors, in order to rename them and treat them as populations of the same species. Multiple combinations were found for four species (V. alutacea, V. epidendroides, C. desvauxii and C. rotundifolia)

# Now we reload the manually edited dataset.

dispersal_traits <- read.csv(str_c(git_find(),
                                   "/Analyses/Trait_Variation/checked_dispersal_traits.csv"))

# Since we have three different dispersal syndromes in the database (i.e., anemochory, autochory and zoochory), we are assessing the probability of each dispersal syndrome separately. Thus, we are creating new columns indicating whether a given species presents (Coded as 1) or not (coded as 0) each dispersal syndrome.

dispersal_traits <- dispersal_traits %>%
  mutate(Anemochory = case_when(Dispersal_syndrome == "Anemochory" ~ 1,
                                TRUE ~ 0),
         Autochory = case_when(Dispersal_syndrome == "Autochory" ~ 1,
                               TRUE ~ 0),
         Zoochory = case_when(Dispersal_syndrome == "Zoochory" ~ 1,
                              TRUE ~ 0)) %>%
  column_to_rownames("Species_acceptedLCVP")

# We create the object for pruning the phylogeny.

dispersal_traits_names <- names_for_pruning(dispersal_traits)

# We load our dataset with the populations that we must add for these models...

populations_dispersal_traits <- read.csv(str_c(git_find(),
                                                "/Analyses/Trait_Variation/populations_dispersal_traits.csv"))
populations_dispersal_traits <- as.matrix(populations_dispersal_traits)
colnames(populations_dispersal_traits) <- NULL
populations_dispersal_traits

# (...),add them to the full phylogeny ...

dispersal_traits_tree <- add.taxa.phylo(full_tree, populations_dispersal_traits)

# ... and prune the tree.

dispersal_traits_tree2 <- prune.sample(dispersal_traits_names, phylo = dispersal_traits_tree$tree)
dispersal_traits_tree$node.label <- NULL

# We arrange our dataset and phylogeny to that the species names are in the same order.

cbind(row.names(dispersal_traits), dispersal_traits_tree2$tip.label)

dispersal_traits <- dispersal_traits[order(match(rownames(dispersal_traits), 
                                                 dispersal_traits_tree2$tip.label)), , drop = FALSE]

cbind(row.names(dispersal_traits), dispersal_traits_tree2$tip.label)

# Finally, we run the mode for each dispersal syndrome, starting with anemochory.

anemochory_plr <- phyloglm(Anemochory ~ Growth_form + Distribution + Microhabitat,
                           data = dispersal_traits, phy = dispersal_traits_tree2)

summary(anemochory_plr)

anemochory_plr2 <- phyloglm(Anemochory ~ Growth_form * Distribution * Microhabitat,
                            data = dispersal_traits, phy = dispersal_traits_tree2)

summary(anemochory_plr2)

# Then, with autochory.

autochory_plr <- phyloglm(Autochory ~ Growth_form + Distribution + Microhabitat,
                          data = dispersal_traits, phy = dispersal_traits_tree2)

summary(autochory_plr)

autochory_plr2 <- phyloglm(Autochory ~ Growth_form * Distribution * Microhabitat,
                           data = dispersal_traits, phy = dispersal_traits_tree2)

summary(autochory_plr2)

# and with zoochory

zoochory_plr <- phyloglm(Zoochory ~ Growth_form + Distribution + Microhabitat,
                         data = dispersal_traits, phy = dispersal_traits_tree2)

summary(zoochory_plr)

zoochory_plr2 <- phyloglm(Zoochory ~ Growth_form * Distribution * Microhabitat,
                          data = dispersal_traits, phy = dispersal_traits_tree2)

summary(zoochory_plr2)

#### 9. PLR for dispersal season ####

# To conclude, we are repeating the same routine with seed dispersal season.

# First, we create the seed dispersal season database.

season_traits <- traits_joined %>%
  dplyr::select(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat, Dispersal_period) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Dispersal_period),
         Growth_form != "Liana",
         Growth_form != "Succulent",
         Growth_form != "Tree") %>%
  group_by(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat, Dispersal_period) %>%
  distinct()

write.csv(season_traits, "check_season_traits.csv")

season_traits <- read.csv(str_c(git_find(),
                                "/Analyses/Trait_Variation/checked_season_traits.csv"))

# Now, we are recoding the database, so that we create four columns (one for each
# dispersal season) and record whether the species disperse their seeds during
# each of them.

season_traits <- season_traits %>%
  distinct() %>%
  mutate(ED = case_when(Dispersal_period == "ED" ~ 1,
                        TRUE ~ 0),
         LD = case_when(Dispersal_period == "LD" ~ 1,
                        TRUE ~ 0),
         ER = case_when(Dispersal_period == "ER" ~ 1,
                        TRUE ~ 0),
         LR = case_when(Dispersal_period == "LR" ~ 1,
                        TRUE ~ 0)) %>%
  group_by(Species_acceptedLCVP, Growth_form, Distribution, Microhabitat) %>%
  summarise(ED = max(ED),
            LD = max(LD),
            ER = max(ER),
            LR = max(LR)) %>%
  column_to_rownames("Species_acceptedLCVP")

# We prune the phylogeny.

season_traits_names <-  names_for_pruning(season_traits)

season_traits_tree <- prune.sample(season_traits_names, phylo = full_tree)
season_traits_tree$node.label <- NULL

# Check that the species name are in the same order than in the phylogeny.

cbind(row.names(season_traits), season_traits_tree$tip.label)
season_traits <- season_traits[order(match(rownames(season_traits), 
                                           season_traits_tree$tip.label)), , drop = FALSE]
cbind(row.names(season_traits), season_traits_tree$tip.label)

# Now with these two objects created, we proceed to run the model for each dispersal
# season.

# We start with early-dry season dispersal.

ED_plr <- phyloglm(ED ~ Growth_form + Distribution + Microhabitat,
                   data = season_traits, phy = season_traits_tree)

summary(ED_plr)

ED_plr2 <- phyloglm(ED ~ Growth_form * Distribution * Microhabitat,
                    data = season_traits, phy = season_traits_tree)

summary(ED_plr2)

# Then with late-dry season dispersal.

LD_plr <- phyloglm(LD ~ Growth_form + Distribution + Microhabitat,
                   data = season_traits, phy = season_traits_tree)

summary(LD_plr)

LD_plr2 <- phyloglm(LD ~ Growth_form * Distribution * Microhabitat,
                    data = season_traits, phy = season_traits_tree)

summary(LD_plr2)

# Early-rain dispersal.

ER_plr <- phyloglm(ER ~ Growth_form + Distribution + Microhabitat,
                   data = season_traits, phy = season_traits_tree)

summary(ER_plr)

ER_plr2 <- phyloglm(ER ~ Growth_form * Distribution * Microhabitat,
                    data = season_traits, phy = season_traits_tree)

summary(ER_plr2)

# And finally, late-rain dispersal.

LR_plr <- phyloglm(LR ~ Growth_form + Distribution + Microhabitat,
                   data = season_traits, phy = season_traits_tree)

summary(LR_plr)

LR_plr2 <- phyloglm(LR ~ Growth_form * Distribution * Microhabitat,
                    data = season_traits, phy = season_traits_tree)

summary(LR_plr2)