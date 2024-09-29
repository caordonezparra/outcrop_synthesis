# This code was designed to assess whether seven seed functional traits exhibit a significant phylogenetic signal. In the case of quantitative traits, it also evaluates the depth of such signal, testing whether its significance stands across different taxonomic levels.

# Authors: Carlos A. Ordóñez-Parra (carlos.ordonez.parra@gmail.com) 
#          and Daniel Negreiros (negreiros.eco@gmail.com)
# Last updated in September 28th, 2024.

#### 0. Loading packages (assuming they are already installed) #####

library(gert) # v. 2.0.1 for managing Git repositories.
library(stringr) # v. 1.5.1 for string manipulation.
library(readr) # readr v. 2.1.5. for loading the data sets.
library(dplyr) # dplyr v. 1.1.4 for dataset manipulation.
library(tidyr) # tidyr v. 1.3.0 for cleaning up the dataset.
library(tibble) # tibble v. 3.2.1 for dataframe manipulations.
library(phytools) # phytools v. 2.1-1 for phylogenetic signal testing.
library(picante) # picante v. 1.8.2 for phylogeny manipulation.
library(car) # car v. 3.1-2 for logit transformation.
library(ade4) # ade4 v. 1.7-19 for creating the objects required for Pavoine et al. (2010) approach.
library(adiv) # adiv v. 2.2 for assessing the phylogenetic signal of qualitative traits
library(LCVP) # 1.0.4. for the Leipzig Catalogue of Vascular Plants.
library(lcvplants) # 2.0 for functions to standardize taxonomy following the LCVP.

#### 1. Set working directory ####

setwd(str_c(git_find(),"/Analyses/Phylogenetic_Signal"))

#### 2. Loading datasets ####

# For this code, we will need to datasets. First, we are loading out trait dataset from Rock n' Seeds database (https://doi.org/10.1002/ecy.3852).

traits <- read_csv(str_c(git_find(),"/traits.csv"))

# Second, we load the phylogeny we created (see Analyses/Phylogeny/phylogeny.R).

full_tree <- read.tree(str_c(git_find(),"/Analyses/Phylogeny/outcrop_phylo.tre"))

#### 3. Filtering the trait database ####

# First, we are obtaining the species accepted name following LCVP (Freiberg et al. 2020), given that need to match the names in our trait dataset to the names in the phylogeny we created (which used this system). We are saving this  object as "lcvp_names". This procedure will also allow us to level up subspecies and varieties names to the species level and will aid fitering out taxa identified up to genus. We are also adding manually the name of Pleroma marumbiensis, manually. This species is not present in the LCVP database, but in Flora and Funga do Brasil (http://floradobrasil.jbrj.gov.br/; Last checked, June 16th,  2022). 

# While the accepted name according to Flora do Brasil would be P. marumbiensis, All Tibouchina and Pleroma species where renamed as Pleroma. Therefore, we kept the name of this species as Tibouchina marumbiensis.

lcvp_names <- lcvp_search(traits$Species_reported) %>%
  dplyr::select(Search, Output.Taxon) %>%
  separate(Output.Taxon, into = c("Species", "Genus"), sep = " ", extra = FALSE) %>%
  unite(col = Species_acceptedLCVP, Species, Genus, sep = "_") %>%
  rename(Species_reported = Search) %>%
  distinct() %>%
  mutate(Species_acceptedLCVP = case_when(
    Species_reported == "Pleroma marumbiense" ~ "Tibouchina_marumbiensis",
    TRUE ~ Species_acceptedLCVP))

# Now, we are going to filter the trait database just to keep the species name as reported in the papers and the seven traits we are interested in (i.e., dispersal period, dispersal syndrome, percentage of embryoless seeds, percentage of viable seeds, seed dry mass, seed water content and primary dormancy). We are also keeping information on species Growth form, geographic distribution and microhabitat.

traits_filtered <- traits %>%
  dplyr::select(Order, Family, Genus, Species_reported, Growth_form, Distribution, Microhabitat,
                Dispersal_period, Dispersal_syndrome, Embryoless, Viable, Dry_mass, Water_content,
                Dormancy)

# Now, we join this filtered trait dataset with the table with the accepted names following the LCVP.

traits_joined <- traits_filtered %>%
  left_join(lcvp_names)

# To ensure that the analyses can be reproduced, we have included the traits_joined file that we used in the original analyses as a .csv file. This can be accessed in the Phylogenetic_Signal folder.

traits_joined <- read_csv((str_c(git_find(),"/Analyses/Phylogenetic_Signal/traits_joined.csv")))

#### 4. Phylogenetic signal for seed mass ####

# First, we filter the traits database to keep only the records with seed dry mass values and that have an accepted name.

seed_mass <- traits_joined %>%
  dplyr::select(Order, Family, Genus, Species_acceptedLCVP, Dry_mass) %>%
  filter(!is.na(Species_acceptedLCVP), # Removes Species_accepted name = "NA" (i.e., taxa identified only to genus)
         !is.na(Dry_mass)) %>% # Removes observations with no mass values.
  group_by(Order, Family, Genus, Species_acceptedLCVP) %>% # Since we might have more than one mass
  # value for a single species, we are calculating the mean dry mass.
  summarise(Mass = mean(Dry_mass)) %>%
  column_to_rownames(var = "Species_acceptedLCVP") # Finally, we are putting the species
  # names as the dataset row names to meet the format of the analysis.

# Second, we are extracting our species names and create an objetct that will allow us to prune the phylogenetic tree and keep only the species with dry seed  mass data. Since this is a procedure that will be done repeteadly, we are creating a function for that, called "names_for_pruning".

# The function takes the dataframe (x) and create an object with all species names.Then, it will take these names and use it as headers of the dataset and remove all observations in the dataset. At the end, our dataset will have as many columns as species in the seed mass dataset, and zero rows.

names_for_pruning <- function(x) {
  names_pruning <- as.data.frame(t(rownames(x)))
  names(names_pruning) <- names_pruning[1,]
  names_pruning <- names_pruning[-1,]
}

# Now that we defined the function, we are going to use it to create the object with the names we need to prune the phylogeny.

seed_mass_names <- names_for_pruning(seed_mass)

# Now, we pruning our phylogenetic tree to keep only the species present in the seed mass dataset.

mass_tree <- prune.sample(seed_mass_names, phylo = full_tree)

# Now, for the phylogenetic signal test, we need that the order in which species appear in the trait dataset is the same as in the phylogeny. We can do that using the following line:

seed_mass <- seed_mass[order(match(rownames(seed_mass), mass_tree$tip.label)), , drop = FALSE]

# Before conducting the phylogenetic signal test itself, we are testing whether the data follow or not a normal distribution. Consiering that we have > 200 observations of seed mass, we will test the normality of the data using the Kolmogorov-Smirnov test.

ks.test(seed_mass$Mass, pnorm, mean(seed_mass$Mass), sd(seed_mass$Mass))
hist(seed_mass$Mass)

# Since raw data are not normal, we are log transforming our data to assess whether this improves their normality. We will see that with this transformation, the data follow a normal distribution.

ks.test(log(seed_mass$Mass), pnorm, mean(log(seed_mass$Mass)), sd(log(seed_mass$Mass)))
hist(log(seed_mass$Mass))

# After resolving that, we will run the phylogenetic signal test using the log-transformed data.

mass_lambda <- phylosig(tree = mass_tree, x = log(seed_mass$Mass), 
                  method = "lambda", test = TRUE)
mass_lambda

b# Now that we found that seed mass exhibits a significant phylogenetic signal, we will test until which taxonomic level it remains significant. For that purpose, we will be using Moran's I as implemented in the correlogram.formula of the ape package.

# First, we are creating an object where the species of the seed mass dataset appear as a column (i.e. not as rownames) and we will log-transform mass values.

mass_Moran_data <- rownames_to_column(seed_mass, var = "Species_acceptedLCVP") %>%
  mutate(LogMass = log(Mass))

# Now, we apply the function, in order to calculate Moran's I for the taxonomic levels available in our dataset: Genus, Family and Order

mass_Moran_test <- correlogram.formula(LogMass ~ Order/Family/Genus,
                                       data = mass_Moran_data)

mass_Moran_test
plot(mass_Moran_test)

# Based on this test, we can say that seed mass a significant autocorrelation across all taxonomic levels, being positive for Genus and Family, but negative for Order.

#### 5. Phylogenetic signal for seed water content ####

# Now, we will repeate the same procedure we did with seed mass data, but now for seed water content.

# Create the dataset with the water content value for each species.

water_content <- traits_joined %>%
  dplyr::select(Order, Family, Genus, Species_acceptedLCVP, Water_content) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Water_content)) %>%
  group_by(Order, Family, Genus, Species_acceptedLCVP) %>%
  summarise(Water_content = mean(Water_content/100)) %>% # The data in the database
  # is stored as percentages. We are dividing it by 100 to make it a proportion 
  # for the transformation we are using later (see below).
  column_to_rownames( var = "Species_acceptedLCVP")

# Creating the object with the species names to prune the phylogenetic tree.

water_content_names <- names_for_pruning(water_content)

# Pruning the phylogenetic tree.

water_tree <- prune.sample(water_content_names, phylo = full_tree)

# Reordering the water content database to match the order of the phylogenetic tree

water_content <- water_content[order(match(rownames(water_content), water_tree$tip.label)), , drop = FALSE]

# Test the normality of the water content data.

ks.test(water_content$Water_content, pnorm, mean(water_content$Water_content), sd(water_content$Water_content))
hist(water_content$Water_content)

# Since raw data did not followed a normal distribution, and water content is a a proportion, we are using an angular transformation to try to fit the data into a normal distribution.

ks.test(asin(sqrt(water_content$Water_content)), pnorm, mean(asin(sqrt(water_content$Water_content))), sd(asin(sqrt(water_content$Water_content))))
hist(asin(sqrt(water_content$Water_content)))

# The angular transformation improved the normality of the data. Yet, they still (although marginally) from a normal distribution. Therefore, we are testing a  logit transformation.

ks.test(logit(water_content$Water_content), pnorm, mean(logit(water_content$Water_content)), sd(logit(water_content$Water_content)))
hist(logit(water_content$Water_content))

# Considering the results of our normality test, we are running the phylogenetic signal test with the logit transformed data.

water_test <- phylosig(tree = water_tree, x = logit(water_content$Water_content), 
                  method = "lambda", test = TRUE)
water_test

# Now that we found that this trait has a significant phylogenetic signal, we will test until which taxonomic level this signal remains. First, we are creating an object where the required object with the water content dataset.

water_Moran_data <- rownames_to_column(water_content, var = "Species_acceptedLCVP") %>%
  mutate(LogitWater = logit(Water_content))

# Now, we apply the function, in order to calculate Moran's I for the taxonomic levels available in our dataset: Genus, Family and Order

water_Moran_test <- correlogram.formula(LogitWater ~ Order/Family/Genus,
                                       data = water_Moran_data)

water_Moran_test
plot(water_Moran_test)

# The test indicates that water content exhibits a significant and positive autocorrelation at the Genus and Family level. Still, no significant autocorrelation for the Order level was detected.

#### 6. Phylogenetic signal for the percentage of embryoless seeds ####

# Now, we will proceed with the same routine for the embryoless seed data. First, we'll create the dataset with the mean proportion of embryoless seeds for each species.

embryoless <- traits_joined %>%
  dplyr::select(Order, Family, Genus, Species_acceptedLCVP, Embryoless) %>%
  filter(!is.na(Species_acceptedLCVP),
         Species_acceptedLCVP != "NA_NA",
         !is.na(Embryoless)) %>%
  group_by(Order, Family, Genus, Species_acceptedLCVP) %>%
  summarise(Embryoless = mean(Embryoless/100)) %>% # As with water content, we
  # are going to save this variable as proportion, rather than a percentage.
  column_to_rownames(var = "Species_acceptedLCVP")

# Extracting the species names.

embryoless_names <- names_for_pruning(embryoless)

# Pruning the phylogenetic tree to keep only those species with data on embryoless seed percentage.

embryoless_tree <- prune.sample(embryoless_names, phylo = full_tree)

# Reordering the embryoless seed dataset so it has the same order of the phylogenetic tree.

embryoless <- embryoless[order(match(rownames(embryoless), embryoless_tree$tip.label)), , drop = FALSE]

# Testing the normality of the embryoless seed data. 

ks.test(embryoless$Embryoless, pnorm, mean(embryoless$Embryoless), sd(embryoless$Embryoless))
hist(embryoless$Embryoless)

# Since raw data did not follow a normal distribution, we are using an angular transformation

ks.test(asin(sqrt(embryoless$Embryoless)), pnorm, mean(asin(sqrt(embryoless$Embryoless))), sd(asin(sqrt(embryoless$Embryoless))))
hist(asin(sqrt(embryoless$Embryoless)))

# The angular transformation did improve the normality of the data. Yet, the distribution of the data still differs from a normal distribution. So, we are testing a logit transformation.

ks.test(logit(embryoless$Embryoless), pnorm, mean(logit(embryoless$Embryoless)), sd(logit(embryoless$Embryoless)))
hist(logit(embryoless$Embryoless))

# With the logit transformation the distribution is only marginally different from a normal distribution. Therefore, we are using it for our phylogenetic signal test.

embryoless_test <- phylosig(tree = embryoless_tree, x = logit(embryoless$Embryoless), 
                  method = "lambda", test = TRUE)
embryoless_test

# Now that we found that this trait has a significant phylogenetic signal, we will test until which taxonomic level this signal remains.

# First, we are creating an object where the required object with the embryoless percentage dataset

embryoless_Moran_data <- rownames_to_column(embryoless, var = "Species_acceptedLCVP") %>%
  mutate(LogitEmbryoless = logit(Embryoless))

# Now, we apply the function, in order to calculate Moran's I for the taxonomic levels available in our dataset: Genus, Family and Order

embryoless_Moran_test <- correlogram.formula(LogitEmbryoless ~ Order/Family/Genus,
                                        data = embryoless_Moran_data)

embryoless_Moran_test
plot(embryoless_Moran_test)

# The dataset indicates that the percentage of embryoless seeds exhibits a significant, positive phylogenetic signal at the Genus and Family level. No significant autocorrelation was detected at Order level.

#### 7. Phylogenetic signal for the percentage of viable seeds ####

# To conclude with the phylogenetic signal tests for our quantitative traits, we are repeating our routine with the data on the percentage of viable seeds.

# We filter out the trait database to keep the species with data for this trait.

viable <- traits_joined %>%
  dplyr::select(Order, Family, Genus, Species_acceptedLCVP, Viable) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Viable)) %>%
  group_by(Order, Family, Genus, Species_acceptedLCVP) %>%
  summarise(Viable = mean(Viable/100)) %>% # We are transforming the percentage into a proportion.
  column_to_rownames(var = "Species_acceptedLCVP")

# We are creating the objetct with the species names present in the previous dataset.

viable_names <- names_for_pruning(viable)

# Now, we prune our phylogenetic tree to keep only the species with data on seed viability.

viable_tree <- prune.sample(viable_names, phylo = full_tree)

# Then, we arrange the dataset so that species appear in the same order as in the phylogeny.

viable <- viable[order(match(rownames(viable), viable_tree$tip.label)), , drop = FALSE]

# Now, we test if the raw viability data follow a normal distribution or not.

ks.test(viable$Viable, pnorm, mean(viable$Viable), sd(viable$Viable))
hist(viable$Viable)

# Since raw data do not follow a normal distribution, we are trying an angular transformation.

ks.test(asin(sqrt(viable$Viable)), pnorm, mean(asin(sqrt(viable$Viable))), sd(asin(sqrt(viable$Viable))))
hist(asin(sqrt(viable$Viable)))

# The angular transformation improved normality, but yet the distribution still differs. So, we are now trying with a logit transformation.

ks.test(logit(viable$Viable), pnorm, mean(logit(viable$Viable)), sd(logit(viable$Viable)))
hist(logit(viable$Viable))

# Since logit-transformed data follow a normal distribution, we are using this transformation for the estimation of the phylogenetic singal.

viable_test <- phylosig(tree = viable_tree, x = logit(viable$Viable), 
                  method = "lambda", test = TRUE)
viable_test

# Now that we found that this trait has a significant phylogenetic signal, we will test until which taxonomic level this signal remains.

# First, we are creating an object where the required object with the water content dataset

viable_Moran_data <- rownames_to_column(viable, var = "Species_acceptedLCVP") %>%
  mutate(LogitViable = logit(Viable))

# Now, we apply the function, in order to calculate Moran's I for the taxonomic levels available in our dataset: Genus, Family and Order

viable_Moran_test <- correlogram.formula(LogitViable ~ Order/Family/Genus,
                                             data = viable_Moran_data)

viable_Moran_test
plot(viable_Moran_test)

# The test indicates that the percentage of viable seeds exhibit a significant and positive autocorrelation at the Genus and Family level, but not at the Order one.

#### 8. Phylogenetic signal for dispersal syndrome ####

# For the phylogenetic signal tests for categorical traits we are using the the approach of Pavoine et al. (2010) as implemented in the adiv package.

# The first trait we are testing is the dispersal syndrome, a categorial trait with three levels: anemochory, zoochory and autochory. To start, we create a database with the species and the dispersal syndrome. For this trait, each species could only be assigned to a single state.

syndrome <- traits_joined %>%
  dplyr::select(Species_acceptedLCVP, Dispersal_syndrome) %>%
  filter(!is.na(Species_acceptedLCVP), # Filtering out species identified until Genus
         !is.na(Dispersal_syndrome)) %>% # and with no record of dispersal syndrome.
  distinct() %>% # Getting only unique combinations.
  column_to_rownames(var = "Species_acceptedLCVP")

# Now, we create the object with the species names present in this dataset.

syndrome_names <- names_for_pruning(syndrome)

# Then, we prune our phylogenetic tree to keep only the species with data on dispersal syndrome.

syndrome_tree <- prune.sample(syndrome_names, phylo = full_tree)

# Now, we are rearranging the dataset, so that species appear in the same order as the phylogeny.

syndrome <- syndrome[order(match(rownames(syndrome), syndrome_tree$tip.label)), , drop = FALSE]

# Now that we are sure they are in the same order, we proceed to creat the objetcs we need for the test.

# To start, we generate our distance matrix. The first step is to create a ktab object, using the ktab.list.df function from the ade4 package.

ktab_syndrome <- ktab.list.df(list(syndrome))

# Then, we use the dist.ktab function to create the distance matrix. Since this is a categorical trait with more than two possible states, we specificy this in 'type' argument using "N" (i.e., nominal).

dist_syndrome <- dist.ktab(ktab_syndrome, type = "N")

comm.Syndrome <- setNames(rep(1, nrow(syndrome)), row.names(syndrome))

# With this two objects, we proceed to perfoming the analysis.

syndrome_test <- rtestdecdiv(phyl = syndrome_tree, comm.Syndrome, 
                              dist_syndrome, nrep = 9999,
                              ties.method = "average", vranking = "droot")

syndrome_test

# The first test indicates that dispersal syndrome exhibits a significant root-to-tip skewness, implying the presence of a significant phylogenetic signal.

#### 9. Phylogenetic signal for dormancy ####

# Now we proceed to repeat the same routine with the dormancy trait. This trait can take two values: ND (non-dormant) and D (dormant). Yet, different populations  of a given species can have different values for this trait. Therefore, we are treating this trait as a multichoice variable.

dormancy <- traits_joined %>%
  dplyr::select(Species_acceptedLCVP, Dormancy) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Dormancy),
         Dormancy != "NC") %>% # Excluding species observations with NC 
  # (i.e., non-conclusive) records
  distinct() %>%
  mutate(Record = 1) %>% # We are creating a new column named "Record" whose value 
  # is always 1. This will indicate the presence of a given state in a species.
  pivot_wider(names_from = Dormancy, values_from = Record, values_fill = 0) %>%
  # Now we will divide the Dormancy into two columns, one with ND and the other
  # with ND and keeping the values (i.e., the 1s) from the Record column. This will
  # indicate that species has a record for a given dormancy trait state. When
  # there are no values available, we are filling the blanks with 0, to indicate
  # the absence of records.
  column_to_rownames(var = "Species_acceptedLCVP")

# From this point, the routine is overall the same as with other traits, with just one exception that is going to be pointed out below.

# Create the object with the species names present in the dormancy dataset.

dormancy_names <- names_for_pruning(dormancy)

# Prune the phylogenetic tree to keep only the species with data on seed dormancy.

dormancy_tree <- prune.sample(dormancy_names, phylo = full_tree)

# Rearrange the names in the trait database to match the order they appear in the phylogeny.

dormancy <- dormancy[order(match(rownames(dormancy), dormancy_tree$tip.label)), , drop = FALSE]

# Now we build the distance matrix. The exceptions that were previously mentioned are here. 

# Since we have a multichoice variable we need to perform a previous step before creating the ktab objetc. We will use the prep.binary to create an object that will keep the information from our different columns and take them to the distance calculation. In this case, we set the argument col.blocks to 2, since we have two possibilities (i.e., ND or D).

prep_dormancy <- prep.binary(dormancy, col.blocks = 2)

# Now we proceed to create the ktab object.

ktab_dormancy <- ktab.list.df(list(prep_dormancy))

# And create the distance matrix. Please note that here type is set to "B", which stands for "binomial".

dist_dormancy <- dist.ktab(ktab_dormancy, type = "B")
comm.Dormancy <- setNames(rep(1, nrow(dormancy)), row.names(dormancy))

# Now that we have created these two objects, we proceed to run the trait diversity decomposition test.

dormancy_test <- rtestdecdiv(dormancy_tree, comm.Dormancy, 
                              dist_dormancy, nrep = 9999,
                              ties.method = "average", vranking = "droot")

dormancy_test

# Based on the result of the Test 1, we can say that dormancy exhibits a significant root-to-tip skewness, indicating the presence of a significant phylogenetic signal.

#### 10. Phylogenetic signal for dispersal season ####

# To conclude our phylogenetic signal tests, we proceed to run the trait diversity decomposition test with the last trait: dispersal season.

# This is also a categorical trait, with four posible states: Early-rain (ER), late-rain (LR), early-dry (ED) and late-dry (LD) season. Yet, as happens with dormancy, species can have records of dispersal in more than one season. Therefore, we will also treating this as a multichoice variable.

# We will preparing the dataset using the same procedure as with dormancy.

season <- traits_joined %>%
  dplyr::select(Species_acceptedLCVP, Dispersal_period) %>%
  filter(!is.na(Species_acceptedLCVP),
         !is.na(Dispersal_period)) %>%
  distinct() %>%
  mutate(Record = 1) %>%
  pivot_wider(names_from = Dispersal_period, values_from = Record, values_fill = 0) %>% 
  column_to_rownames(var = "Species_acceptedLCVP")

# We are creating the object with the species names present in the previous dataset.

season_names <- names_for_pruning(season)

# We prune our phylogenetic tree to keep only the species with data on dispersal season.

season_tree <- prune.sample(season_names, phylo = full_tree)

# Rearranging the trait dabase so the order of species is the same as the phylogenetic tree.

cbind(row.names(season), season_tree$tip.label)
season <- season[order(match(rownames(season), season_tree$tip.label)), , drop = FALSE]
cbind(row.names(season), season_tree$tip.label)

# We generate the distance matrix. Please note that for the preparing archive, we should set the col.blocks to 4 (i.e., considering the four possible states for the trait).

prep_season <- prep.binary(season, col.blocks = 4)
ktab_season <- ktab.list.df(list(prep_season))
dist_season <- dist.ktab(ktab_season, type = "B")

comm.Season <- setNames(rep(1, nrow(season)), row.names(season))

# With these two objects created, we proceed to run the trait diversity decomposition test.

season_test <- rtestdecdiv(season_tree, comm.Season, 
                              dist_season, nrep = 9999,
                              ties.method = "average", vranking = "droot")

season_test

# According to the test, dispersal season exhibits a significant root-to-tip skewness, implying a significant phylogenetic signal.