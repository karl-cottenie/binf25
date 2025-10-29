#Yumna Khan

###***************************
## MBINF 6210
##
## Yumna Khan - 1094825
##
## 2025-10-23
##
##***************************

## _ Packages used -------
library(rentrez)
library(Biostrings)
library(tidyverse)
library(dplyr)
library(randomForest)


## _ Search and Read Data -------

df_anopheles_search <- entrez_search(db = "nuccore", term = "anopheles [ORGN] AND COI [gene]", use_history = T)
df_anopheles_search

# The search returns 10730 hits! On NCBI, in nucleotides data base, search up the term and download the fasta file

# Read the file
st_anopheles <- readDNAStringSet("../data/sequence.fasta")
st_anopheles

# Create the dataframe with sequence titles and DNA strings
df_anopheles <- data.frame(title = names(st_anopheles), sequence = paste(st_anopheles))
head(df_anopheles)

# TODO: Explore data more (i.e. make histogram)


## _ Obtain Nucleotides -------

# Create a new column to convert the vector string into a manipulative DNAstring set
df_anopheles$nucleotides2 <- DNAStringSet(df_anopheles$sequence)
class(df_anopheles$nucleotides2)
head(df_anopheles)

# Compute the dinucleotides frequencies and bind as new columns
df_anopheles <- cbind(
  df_anopheles,
  as.data.frame(dinucleotideFrequency(df_anopheles$nucleotides2, as.prob = TRUE))
)
head(df_anopheles)


## _ Add Species and Unique ID to DF -------

# Add the unique id and species name to df
# The species name is the second and third word in the title, while the unique id is the first word in the title
df_anopheles <- df_anopheles %>%
  mutate(species_name = word(title, 2L, 3L)) %>%
  mutate(unique_id = word(title, 1L))
df_anopheles

# Count the unique id
df_anopheles$unique_id %>%
  unique() %>%
  length()

# Remove the nucleotides2 column since it won't be used later
df_anopheles <- df_anopheles %>%
  select(-nucleotides2) %>%
  as_tibble()


## _ Filter Data -------

# Obtain species name in list format if there's 200+ of them
ls_vectors <- df_anopheles %>%
  group_by(species_name) %>%
  count() %>%
  arrange(desc(n)) %>% # View()
  filter(n > 200) %>%
  pull(species_name)
ls_vectors

# List selected species
df_anopheles_vector <-
  df_anopheles %>%
  filter(species_name %in% ls_vectors)

# Check number of each species is 200+
df_anopheles_vector %>%
  group_by(species_name) %>%
  count()

# Remove unsued objects to free memory
rm(df_anopheles, st_anopheles)

# Obtain smallest number of species count
smaller_sample <- min(table(df_anopheles_vector$species_name))
smaller_sample


## _ Get Validation Set -------

# Set random seed for reproducibility
set.seed(51)

# Sample 20% of each species for validation
df_anopheles_validation <- df_anopheles_vector %>%
  group_by(species_name) %>%
  sample_n(floor(0.2 * smaller_sample))
df_anopheles_validation

# Verify that 20% of each species are for validation set
df_anopheles_validation %>%
  group_by(species_name) %>%
  count()

# TODO: Create a test set for final model evaluation


## _ Train Data -------

# Sample 80% of training data, excluding species and unique id from validation set
set.seed(40)
df_anopheles_training <- df_anopheles_vector %>%
  filter(!species_name %in% df_anopheles_validation$unique_id) %>%
  group_by(species_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

names(df_anopheles_training)

# Verify number of training samples per species
df_anopheles_training %>%
  group_by(species_name) %>%
  count()


## _ Random Forest -------

# Use the dinucleotide frequencies to make a RandomForest Model of 100 trees
anopheles_class <- randomForest(
  x = df_anopheles_training[, 3:18],
  y = as.factor(df_anopheles_training$species_name),
  ntree = 100, importance = TRUE
)


## _ Plot -------

# Output of RandomForest Model:
# 1. Mean Decrease Accuracy (decrease accuracy when values are randomly shuffled)
# 2. Mean Decrease Gini (measure of class mixing)
# 3. Error rate as number of trees increase
anopheles_class
plot(anopheles_class)
varImpPlot(anopheles_class)


## _ Validate -------

# Predict species label for validation set using training model
predict_validation <- predict(
  anopheles_class,
  df_anopheles_validation[, c(19, 3:18)]
)
predict_validation

## _ Performance Metric -------

# Table of observed vs predicted (confusion matrix)
table(
  observed = df_anopheles_validation$species_name,
  predicted = predict_validation
)

# TODO:
# - Using test set, test the RandomForest model
# - Compute performance metrics (e.g. F1 score, confusion matrix, mean square error)
# - Visualize model accuracy and feature importance

